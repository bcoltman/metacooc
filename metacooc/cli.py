#!/usr/bin/env python3

## placeholder for now

import argparse
from aastk.pasr import extract_matching_sequences

def parse_cli():
    parser = argparse.ArgumentParser(description="Amino Acid Sequence Toolkit")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for 'pasr'
    pasr_parser = subparsers.add_parser("pasr", help="PASR: protein alignment score ratio")
    pasr_parser.add_argument("blast_tab", help="BLAST/DIAMOND tabular output file")
    pasr_parser.add_argument("seq_file", help="FASTA or FASTQ file")
    pasr_parser.add_argument("out_fasta", help="Output file for matched sequences")

    args = parser.parse_args()

    if args.command == "pasr":
        extract_matching_sequences(args.blast_tab, args.seq_file, args.out_fasta)
