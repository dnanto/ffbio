#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO

from itertools import chain


def qualifiers(record):
    feats = getattr(record, "features", [{}])
    return record, getattr(feats[0], "qualifiers", {})


def parse_argv(argv):
    parser = ArgumentParser(
        description="generate a table of qualifier metadata",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=FileType(), help="the sequence file")
    parser.add_argument("keys", nargs="*")
    parser.add_argument("-default", default="")
    parser.add_argument("-joiner", dest="joi", default=";")
    parser.add_argument("-separator", dest="sep", default="\t", nargs="*")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])
    with args.file as file:
        rec_quals = map(qualifiers, SeqIO.parse(file, "genbank"))

        if not args.keys:
            args.keys = set(chain.from_iterable(item[-1].keys() for item in rec_quals))
            rec_quals = map(qualifiers, SeqIO.parse(file, "genbank"))

        print("accver", *args.keys, sep=args.sep)
        for rec, qual in rec_quals:
            print(
                rec.id,
                *(args.joi.join(qual.get(key, args.default)) for key in args.keys),
                sep=args.sep
            )

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
