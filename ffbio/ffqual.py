#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO

from itertools import chain


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
    parser.add_argument("-default", default="?", help="the default value for missing entries")
    parser.add_argument("-joi", default=";", help="the field value join character")
    parser.add_argument("-sep", default="\t", nargs="*", help="the table separator")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])
    with args.file as stream:
        records = SeqIO.parse(stream, "genbank")

        if not args.qual:
            records = list(records)
            args.qual = set(chain.from_iterable(ele.keys() for ele in map(qualifiers, records)))

        print("accver", *args.qual, sep=args.sep)
        for rec in records:
            obj = qualifiers(rec)
            row = (args.joi.join(obj.get(key, args.default)) for key in args.qual)
            print(rec.id, *row, sep=args.sep)


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
