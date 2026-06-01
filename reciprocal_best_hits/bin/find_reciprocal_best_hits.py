#!/usr/bin/env python3
"""
Find reciprocal best hits (RBH) from forward and reverse BLAST results.

For each query-subject pair, identifies cases where A's best hit is B and
B's best hit is A. Handles ties at the top bitscore gracefully (e.g. from
recent gene duplicates).

Output includes the forward BLAST columns plus reverse alignment stats
(rev_pident, rev_length, rev_gapopen, rev_bitscore) for comparison.
"""

import argparse
import csv
import sys


ALIGN_COLS = ['pident', 'length', 'gapopen', 'bitscore']


def find_col_indices(cols, names):
    """Return a dict of column name -> index for the requested column names."""
    indices = {}
    for name in names:
        if name in cols:
            indices[name] = cols.index(name)
    return indices


def parse_best_hits(filepath, bitscore_idx):
    """Return dict: query -> {score, hits} for all subjects tied at the best bitscore."""
    best = {}
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            qid, sid = row[0], row[1]
            score = float(row[bitscore_idx])

            if qid not in best:
                best[qid] = {'score': score, 'hits': [(sid, row)]}
            elif score > best[qid]['score']:
                best[qid] = {'score': score, 'hits': [(sid, row)]}
            elif score == best[qid]['score']:
                best[qid]['hits'].append((sid, row))
    return best


def find_rbh(forward_file, reverse_file, outfmt_cols, output_file):
    """Identify reciprocal best hits and write results."""
    cols = outfmt_cols.split()
    col_idx = find_col_indices(cols, ALIGN_COLS)

    missing = [c for c in ALIGN_COLS if c not in col_idx]
    if missing:
        sys.exit(
            f"ERROR: --outfmt_cols must include {', '.join(ALIGN_COLS)} "
            f"(missing: {', '.join(missing)})"
        )

    bitscore_idx = col_idx['bitscore']

    forward = parse_best_hits(forward_file, bitscore_idx)
    reverse = parse_best_hits(reverse_file, bitscore_idx)

    # Build header: forward columns + reverse alignment stats
    rev_align_header = ['rev_' + c for c in ALIGN_COLS]
    header = '\t'.join(cols + rev_align_header)

    with open(output_file, 'w') as out:
        out.write(header + '\n')
        for qid in sorted(forward):
            for sid, fwd_row in forward[qid]['hits']:
                if sid in reverse:
                    for rev_sid, rev_row in reverse[sid]['hits']:
                        if rev_sid == qid:
                            rev_vals = [rev_row[col_idx[c]] for c in ALIGN_COLS]
                            out.write('\t'.join(fwd_row + rev_vals) + '\n')
                            break


def main():
    parser = argparse.ArgumentParser(
        description='Find reciprocal best hits from forward and reverse BLAST results.'
    )
    parser.add_argument(
        '--forward', required=True,
        help='Forward BLAST results (query vs subject) in tabular format.'
    )
    parser.add_argument(
        '--reverse', required=True,
        help='Reverse BLAST results (subject vs query) in tabular format.'
    )
    parser.add_argument(
        '--outfmt_cols', required=True,
        help='Space-separated BLAST outfmt column names (must include pident, length, gapopen, bitscore).'
    )
    parser.add_argument(
        '--output', default='reciprocal_best_hits.tsv',
        help='Output file path (default: reciprocal_best_hits.tsv).'
    )
    args = parser.parse_args()

    find_rbh(args.forward, args.reverse, args.outfmt_cols, args.output)


if __name__ == '__main__':
    main()
