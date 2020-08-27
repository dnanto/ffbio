#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import OrderedDict
from itertools import groupby
from pathlib import Path
from signal import SIG_DFL, SIGPIPE, signal

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

from ffbio.ffparse import ffparse, ffsniff


def keyfunc(record):
    return seguid(record.seq)


def parse_argv(argv):
    parser = ArgumentParser(
        description="compute the unique set of sequences",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=FileType(), help="the sequence file")
    parser.add_argument("-map", type=Path, help="the path prefix for the output files")
    parser.add_argument("-fmt", default="fasta", help="the sequence file format")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])

    # aggregate sequences by unique hash
    with args.file as stream:
        keyfunc = lambda val: seguid(val.seq)
        records = OrderedDict(
            (key, list(val))
            for key, val in groupby(sorted(SeqIO.parse(stream, args.fmt), key=keyfunc), keyfunc)
        )

    # map unique sequence index to unique hash
    if args.map:
        with args.map.open("w") as stream:
            width = len(str(len(records)))
            print("idx", "key", "id", "description", "length", sep="\t", file=stream)
            for idx, ele in enumerate(records.items(), start=1):
                idx, key, val = f"{idx:0{width}d}", ele[0], ele[1]
                for rec in val:
                    print(idx, key, rec.id, rec.description, len(rec), sep="\t", file=stream)
                    rec.id = idx

    # output unique sequences
    SeqIO.write((val[0] for val in records.values()), sys.stdout, args.fmt)

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
