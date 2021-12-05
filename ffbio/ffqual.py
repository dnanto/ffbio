#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from itertools import chain
from signal import SIG_DFL, SIGPIPE, signal

from Bio import SeqIO


def qualifiers(rec):
    feat = getattr(rec, "features", [{}])
    return getattr(feat[0], "qualifiers", {})


def parse_argv(argv):
    parser = ArgumentParser(
        description="generate a table of qualifier metadata",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=FileType(), help="the sequence file")
    parser.add_argument("qual", nargs="*", help="the qualifier keys")
    parser.add_argument("-description", action="store_true", help="the flag to include description")
    parser.add_argument("-default", default="?", help="the default value for missing entries")
    parser.add_argument("-joi", default=";", help="the field value join character")
    parser.add_argument("-sep", default="\t", nargs="*", help="the table separator")
    parser.add_argument("-sort", action="store_true", help="the flag to sort by id")
    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])
    with args.file as stream:
        records = SeqIO.parse(stream, "genbank")

        if not args.qual:
            records = list(records)
            args.qual = set(chain.from_iterable(ele.keys() for ele in map(qualifiers, records)))

        hdr = ["id"] + (["description"] if args.description else [])
        records = sorted(records, key=lambda record: record.id) if args.sort else records
        print(*hdr, *args.qual, sep=args.sep)
        for rec in records:
            obj = qualifiers(rec)
            row = (args.joi.join(obj.get(key, args.default)) for key in args.qual)
            print(*(getattr(rec, name) for name in hdr), *row, sep=args.sep)


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
