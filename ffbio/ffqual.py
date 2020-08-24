#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO

from itertools import chain


def qualifiers(rec):
    feat = getattr(rec, "features", [{}])
    return getattr(feat[0], "qualifiers", {})


def process_record(rec, args):
    annos = ()
    quals = ()

    if args.anno:
        anno = rec.annotations
        annos = (anno.get(key, args.default) for key in args.anno)
        annos = (args.joi.join(map(repr, val)) if isinstance(val, list) else val for val in annos)

    if args.qual:
        qual = qualifiers(rec)
        quals = (args.joi.join(qual.get(key, args.default)) for key in args.qual)

    return (rec.id, *annos, *quals)


def parse_argv(argv):
    parser = ArgumentParser(
        description="generate a table of qualifier metadata",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=FileType(), help="the sequence file")
    parser.add_argument("-anno", nargs="+", default=[], help="the annotation keys")
    parser.add_argument("-qual", nargs="+", default=[], help="the qualifier keys")
    parser.add_argument(
        "-all-anno", action="store_true", default=[], help="the flag to get all annotation keys"
    )
    parser.add_argument(
        "-all-qual", action="store_true", default=[], help="the flag to get all qualifier keys",
    )
    parser.add_argument("-default", default="?", help="the default value for missing entries")
    parser.add_argument("-joiner", dest="joi", default=";", help="the field value join character")
    parser.add_argument(
        "-separator", dest="sep", default="\t", nargs="*", help="the table separator"
    )

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])
    with args.file as file:
        records = list(SeqIO.parse(file, "genbank"))

    if args.all_anno:
        args.anno = set(chain.from_iterable(ele.annotations.keys() for ele in records))
    if args.all_qual:
        args.qual = set(chain.from_iterable(ele.keys() for ele in map(qualifiers, records)))

    print("accver", *args.anno, *args.qual, sep=args.sep)
    for rec in records:
        print(*process_record(rec, args), sep=args.sep)

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
