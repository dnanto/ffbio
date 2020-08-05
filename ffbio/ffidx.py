#!/usr/bin/env python3

import os
import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL
from Bio import SeqIO


def parse_argv(argv):
    parser = ArgumentParser(
        description="retrieve records from an indexed set of sequence files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("path", type=Path, help="the sequence flat-file or index path")
    parser.add_argument("-filenames", nargs="+", help="the list of sequence files to index")
    parser.add_argument("-dump", action="store_true", help="the flag to dump all of the records")
    parser.add_argument(
        "-descriptions",
        "-headers",
        action="store_true",
        help="the flag to only output the descriptions",
    )
    parser.add_argument("-entry", nargs="+", help="the accessions to retrieve")
    parser.add_argument("-entry-batch", help="the file of accessions to retrieve")
    parser.add_argument("-fi", default="fasta", help="the sequence file format (input)")
    parser.add_argument("-fo", default="fasta", help="the sequence file format (output)")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])

    if args.path.name.endswith(".db"):
        args.fi = None if args.path.exists() else args.fi
        db = SeqIO.index_db(str(args.path))
    else:
        args.path = sys.stdin if args.path.name == "-" else args.path
        db = SeqIO.to_dict(SeqIO.parse(args.path, args.fi))

    keys = []
    if args.dump:
        keys = db.keys()
    else:
        if args.entry:
            keys = args.entry
        if args.entry_batch:
            with args.entry_batch as stream:
                keys += list(map(str.strip, stream))

    records = (db[key] for key in keys)
    if args.descriptions:
        print(*(record.description for record in records), sep="\n")
    else:
        SeqIO.write(records, sys.stdout, args.fo)

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
