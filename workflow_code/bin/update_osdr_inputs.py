#!/usr/bin/env python
"""Resolve intermediate Amplicon inputs from the OSDR file listing."""

import argparse
import csv
import json
from urllib.parse import urljoin
from urllib.request import urlopen


def file_mapping(data, accession):
    if accession in data and 'files' in data[accession]:
        return {name: item['URL'] for name, item in data[accession]['files'].items()}
    files = data.get('studies', {}).get(accession, {}).get('study_files', [])
    mapping = {}
    for item in files:
        name = item['file_name']
        url = urljoin('https://osdr.nasa.gov', item['remote_url'])
        if name in mapping and mapping[name] != url:
            raise ValueError(f'Multiple URLs for {name}')
        mapping[name] = url
    if not mapping:
        raise ValueError(f'No files found for {accession}')
    return mapping


def resolve(mapping, ending):
    matches = [url for name, url in mapping.items() if name == ending or name.endswith('_' + ending)]
    if len(matches) != 1:
        raise ValueError(f'Expected one OSDR file ending in {ending}; found {len(matches)}. Supply explicit paths in a runsheet.')
    return matches[0]


def update(rows, fields, mapping, entry_point, region):
    suffix = f'_{region}_GLAmpSeq'
    if entry_point == 'asv_tables':
        columns = ['counts_table_path', 'taxonomy_table_path']
        values = [resolve(mapping, f'counts{suffix}.tsv'), resolve(mapping, f'taxonomy{suffix}.tsv')]
        for row in rows:
            row.update(zip(columns, values))
    else:
        columns = ['trimmed_read1_path', 'trimmed_read2_path']
        for row in rows:
            sample = row.get('sample_id') or row.get('Sample Name') or row.get('sample_name')
            paired = (row.get('paired') or row.get('paired_end') or '').lower() == 'true'
            row[columns[0]] = resolve(mapping, f'{sample}{suffix}_R1_trimmed.fastq.gz' if paired else f'{sample}{suffix}_trimmed.fastq.gz')
            row[columns[1]] = resolve(mapping, f'{sample}{suffix}_R2_trimmed.fastq.gz') if paired else ''
    return fields + [name for name in columns if name not in fields], rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--accession', required=True)
    parser.add_argument('--api-url', required=True)
    parser.add_argument('--entry-point', choices=['trimmed_reads', 'asv_tables'], required=True)
    parser.add_argument('--target-region', required=True)
    parser.add_argument('--limit', type=int, default=0)
    args = parser.parse_args()
    try:
        url = args.api_url.replace('*', args.accession).rstrip('/') + '/files/'
        with urlopen(url, timeout=90) as response:
            data = json.load(response)
        with open(args.input, newline='') as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            if args.limit > 0:
                rows = rows[:args.limit]
            fields, rows = update(rows, reader.fieldnames, file_mapping(data, args.accession), args.entry_point, args.target_region)
        with open(args.output, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
