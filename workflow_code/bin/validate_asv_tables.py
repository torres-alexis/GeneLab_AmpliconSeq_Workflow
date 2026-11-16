#!/usr/bin/env python
"""Validate and subset dataset-wide ASV counts and taxonomy to a runsheet."""

import argparse
import csv
from decimal import Decimal, InvalidOperation


def table(path):
    with open(path, newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        rows = list(reader)
    if len(header) < 2 or len(header) != len(set(header)) or any(len(row) != len(header) for row in rows):
        raise ValueError(f'Invalid TSV columns: {path}')
    ids = [row[0] for row in rows]
    if not rows or any(not name for name in ids) or len(ids) != len(set(ids)):
        raise ValueError(f'Empty or duplicate ASV identifiers: {path}')
    return header, rows


def validate(counts, taxonomy, runsheet, out_counts, out_taxonomy):
    with open(runsheet, newline='') as handle:
        samples = [row['sample_id'] for row in csv.DictReader(handle)]
    header, rows = table(counts)
    missing = set(samples) - set(header[1:])
    if missing:
        raise ValueError('Counts table is missing runsheet samples: ' + ', '.join(sorted(missing)))
    indices = [header.index(sample) for sample in samples]
    for row in rows:
        for value in row[1:]:
            try:
                count = Decimal(value)
                valid = count.is_finite() and count >= 0 and count == count.to_integral_value()
            except InvalidOperation:
                valid = False
            if not valid:
                raise ValueError(f'ASV counts must be nonnegative integers: {row[0]}')
    for index in indices:
        if not any(Decimal(row[index]) > 0 for row in rows):
            raise ValueError(f'Sample has no ASV counts: {header[index]}')
    taxa_header, taxa_rows = table(taxonomy)
    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    missing_ranks = [rank for rank in ranks if rank not in taxa_header[1:]]
    if missing_ranks:
        raise ValueError('Taxonomy is missing columns: ' + ', '.join(missing_ranks))
    ranks = ['domain'] + [rank for rank in ('kingdom', 'division') if rank in taxa_header] + ranks[1:]
    taxa_columns = [0] + [taxa_header.index(rank) for rank in ranks]
    taxa_columns += [i for i in range(1, len(taxa_header)) if i not in taxa_columns]
    taxa = {row[0]: row for row in taxa_rows}
    missing = {row[0] for row in rows} - set(taxa)
    if missing:
        raise ValueError('Taxonomy is missing ASVs: ' + ', '.join(sorted(missing)[:10]))
    with open(out_counts, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([header[0]] + samples)
        writer.writerows([row[0]] + [row[i] for i in indices] for row in rows)
    with open(out_taxonomy, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([taxa_header[i] for i in taxa_columns])
        writer.writerows([taxa[row[0]][i] for i in taxa_columns] for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['counts', 'taxonomy', 'runsheet', 'out-counts', 'out-taxonomy']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    try:
        validate(args.counts, args.taxonomy, args.runsheet, args.out_counts, args.out_taxonomy)
    except (ValueError, OSError, StopIteration, KeyError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
