#!/usr/bin/env python3

import os
import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def accverize(index, keys):
    with sqlite3.connect(index) as conn:
        curs = conn.cursor()
        for key in keys:
            curs.execute(
                """
                SELECT key FROM offset_data WHERE key LIKE ? || '.%'
                ORDER BY length(key) DESC, key DESC
                LIMIT 1;
                """,
                (key,),
            )
            result = curs.fetchone()
            yield result[0] if result and len(result) > 0 else key

        curs.close()


def keygetter(db, keys, keyerror=False):
    for key in keys:
        val = db.get(key)
        if val:
            yield val
        elif keyerror:
            raise KeyError


def parse_argv(argv):
    parser = ArgumentParser(
        description="retrieve records from an indexed set of sequence files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("index", help="the index")
    parser.add_argument("-filenames", nargs="+", help="the list of sequence files to index")
    parser.add_argument("-dump", action="store_true", help="the flag to dump all of the records")
    parser.add_argument(
        "-descriptions",
        "-headers",
        action="store_true",
        help="the flag to only output the descriptions",
    )
    parser.add_argument("-entry", nargs="+", help="the accessions to retrieve")
    parser.add_argument(
        "-entry-batch", "--entry-batch", type=FileType(), help="the file of accessions to retrieve"
    )
    parser.add_argument(
        "-keyerror",
        action="store_true",
        help="the flag to exit on key error (if the accession isn't found)",
    )
    parser.add_argument(
        "-xversion",
        action="store_true",
        help="the flag to indicate that the accessions are missing a version",
    )
    parser.add_argument(
        "-fmt-idx", help="the sequence file format of the indexed files, optional if reloading",
    )
    parser.add_argument("-fmt-out", help="the sequence file format (output)")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])

    fmt_idx = args.fmt_idx
    fmt_out = args.fmt_out

    if args.index == "-":
        db = SeqIO.to_dict(SeqIO.parse(sys.stdin, fmt_idx))
    else:
        db = SeqIO.index_db(args.index, filenames=args.filenames, format=fmt_idx)
        if args.index == ":memory:":
            fmt_out = fmt_out or fmt_idx
        else:
            with sqlite3.connect(args.index) as conn:
                meta = dict(conn.execute("SELECT * FROM meta_data"))
                fmt_out = fmt_out or meta["format"]

    keys = []
    if args.dump:
        keys = db.keys()
    else:
        if args.entry:
            keys = args.entry
        if args.entry_batch:
            with args.entry_batch as file:
                keys += list(map(str.strip, file))

    keys = accverize(args.index, keys) if args.xversion else keys
    records = keygetter(db, keys, keyerror=args.keyerror)
    if args.descriptions:
        print(*(record.description for record in records), sep="\n")
    else:
        SeqIO.write(records, sys.stdout, fmt_out)

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
