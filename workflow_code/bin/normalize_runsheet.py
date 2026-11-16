#!/usr/bin/env python
"""Normalize supplied and ISA-derived Amplicon metadata for workflow input."""

import argparse
import csv
import re
from pathlib import Path
from urllib.parse import urlparse


INPUT_COLUMNS = {
    'raw_reads': ('forward', 'reverse'),
    'trimmed_reads': ('trimmed_read1_path', 'trimmed_read2_path'),
    'asv_tables': ('counts_table_path', 'taxonomy_table_path'),
}


def normalize(source, base_dir, sample_column='sample_id', group_column='groups',
              limit=0, force_single_end=None, entry_point='raw_reads'):
    with open(source, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError('Runsheet contains duplicate column names')
        rows = list(reader)
    if not rows:
        raise ValueError('Runsheet contains no samples')
    if entry_point == 'asv_tables' and force_single_end:
        raise ValueError('--force_single_end does not apply to ASV tables')
    if limit > 0:
        rows = rows[:limit]
    ids, read_paths, layouts, tables = set(), set(), set(), set()
    normalized = []
    first, second = INPUT_COLUMNS[entry_point]
    for row in rows:
        if None in row:
            raise ValueError('Runsheet row has more values than column names')
        row = {key: (value or '').strip() for key, value in row.items()}
        sample = row.get(sample_column) if sample_column != 'sample_id' else (
            row.get('sample_id') or row.get('Sample Name') or row.get('sample_name'))
        if not sample or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
            raise ValueError(f'Invalid sample ID: {sample!r}; use letters, numbers, dots, underscores or dashes')
        if sample in ids:
            raise ValueError(f'Duplicate sample ID: {sample}')
        ids.add(sample)
        group = row.get(group_column)
        if not group:
            raise ValueError(f'Missing {group_column} for sample {sample}')
        row.update(sample_id=sample, groups=group)
        if entry_point == 'asv_tables':
            paths = [row.get(first), row.get(second)]
            if not all(paths):
                raise ValueError(f'ASV tables require {first} and {second} for sample {sample}')
        else:
            paired = (row.get('paired') or row.get('paired_end') or '').lower()
            if paired not in ('true', 'false'):
                raise ValueError(f'paired must be true or false for sample {sample}')
            forward = row.get(first) or (row.get('read1_path') if entry_point == 'raw_reads' else None)
            reverse = row.get(second) or (row.get('read2_path') if entry_point == 'raw_reads' else '') or ''
            if not forward or (paired == 'true' and not reverse):
                raise ValueError(f'Missing {first} or {second} for sample {sample}')
            if force_single_end:
                if force_single_end == 'R2':
                    if paired != 'true':
                        raise ValueError(f'Cannot select R2 from single-end sample {sample}')
                    forward = reverse
                paired, reverse = 'false', ''
            row['paired'] = paired
            layouts.add(paired)
            paths = [forward, reverse] if paired == 'true' else [forward]
        resolved = []
        for value in paths:
            location = value if urlparse(value).scheme in ('http', 'https', 'ftp', 's3', 'gs') else str((Path(base_dir) / value).resolve())
            if entry_point != 'asv_tables' and location in read_paths:
                raise ValueError(f'Read file assigned more than once: {location}')
            read_paths.add(location)
            resolved.append(location)
        row[first] = resolved[0]
        row[second] = resolved[1] if len(resolved) == 2 else ''
        if entry_point == 'asv_tables':
            tables.add(tuple(resolved))
        normalized.append(row)
    if len(layouts) > 1:
        raise ValueError('Runsheet must use one read layout per run')
    if len(tables) > 1:
        raise ValueError('Use one count table and one taxonomy table for the full dataset')
    columns = ['sample_id', first]
    if entry_point == 'asv_tables' or 'true' in layouts or second in fields:
        columns.append(second)
    if entry_point != 'asv_tables':
        columns.append('paired')
    columns.append('groups')
    columns += [field for field in fields if field not in columns]
    return columns, normalized


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--base-dir', required=True)
    parser.add_argument('--sample-column', default='sample_id')
    parser.add_argument('--group-column', default='groups')
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--force-single-end', choices=['R1', 'R2'])
    parser.add_argument('--entry-point', choices=INPUT_COLUMNS, default='raw_reads')
    args = parser.parse_args()
    try:
        if Path(args.input).resolve() == Path(args.output).resolve():
            raise ValueError('Output runsheet must be separate from the input runsheet')
        columns, rows = normalize(args.input, args.base_dir, args.sample_column,
                                  args.group_column, args.limit, args.force_single_end, args.entry_point)
        with open(args.output, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
