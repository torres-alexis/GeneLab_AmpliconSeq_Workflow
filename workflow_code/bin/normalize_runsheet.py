#!/usr/bin/env python
"""Normalize supplied and ISA-derived Amplicon metadata for workflow input."""

import argparse
import csv
import re
from pathlib import Path
from urllib.parse import urlparse


def normalize(source, base_dir, sample_column='sample_id', group_column='groups',
              limit=0, force_single_end=None):
    with open(source, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError('Runsheet contains duplicate column names')
        rows = list(reader)
    if not rows:
        raise ValueError('Runsheet contains no samples')
    if limit > 0:
        rows = rows[:limit]
    ids, read_paths, layouts, primers = set(), set(), set(), set()
    normalized = []
    for row in rows:
        if None in row:
            raise ValueError('Runsheet row has more values than column names')
        row = {key: (value or '').strip() for key, value in row.items()}
        sample = row.get(sample_column) or row.get('sample_id') or row.get('Sample Name') or row.get('sample_name')
        if not sample or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
            raise ValueError(f'Invalid sample ID: {sample!r}; use letters, numbers, dots, underscores or dashes')
        if sample in ids:
            raise ValueError(f'Duplicate sample ID: {sample}')
        ids.add(sample)
        group = row.get(group_column)
        if not group:
            raise ValueError(f'Missing {group_column} for sample {sample}')
        paired = (row.get('paired') or row.get('paired_end') or '').lower()
        if paired not in ('true', 'false'):
            raise ValueError(f'paired must be true or false for sample {sample}')
        forward = row.get('forward') or row.get('read1_path')
        reverse = row.get('reverse') or row.get('read2_path') or ''
        if not forward or (paired == 'true' and not reverse):
            raise ValueError(f'Missing read file for sample {sample}')
        if force_single_end:
            if force_single_end == 'R2':
                if paired != 'true':
                    raise ValueError(f'Cannot select R2 from single-end sample {sample}')
                forward = reverse
            paired, reverse = 'false', ''
        layouts.add(paired)
        primers.add((row.get('F_Primer', ''), row.get('R_Primer', '')))
        paths = [forward, reverse] if paired == 'true' else [forward]
        resolved = []
        for value in paths:
            location = value if urlparse(value).scheme in ('http', 'https', 'ftp', 's3', 'gs') else str((Path(base_dir) / value).resolve())
            if location in read_paths:
                raise ValueError(f'Read file assigned more than once: {location}')
            read_paths.add(location)
            resolved.append(location)
        row.update(sample_id=sample, forward=resolved[0], reverse=resolved[1] if paired == 'true' else '', paired=paired, groups=group)
        normalized.append(row)
    if len(layouts) != 1:
        raise ValueError('Runsheet must use one read layout per run')
    if len(primers) != 1:
        raise ValueError('Runsheet contains inconsistent primer sequences')
    columns = ['sample_id', 'forward'] + (['reverse'] if 'true' in layouts or 'reverse' in fields else []) + ['paired', 'groups']
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
    args = parser.parse_args()
    try:
        if Path(args.input).resolve() == Path(args.output).resolve():
            raise ValueError('Output runsheet must be separate from the input runsheet')
        columns, rows = normalize(args.input, args.base_dir, args.sample_column,
                                  args.group_column, args.limit, args.force_single_end)
        with open(args.output, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
