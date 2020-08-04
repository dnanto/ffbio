#!/usr/bin/env python3

import logging
import sqlite3
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from collections import OrderedDict
from datetime import datetime
from itertools import chain
from math import ceil
from os import makedirs
from pathlib import Path
from signal import SIG_DFL, SIGPIPE, signal
from subprocess import PIPE, Popen

from Bio import Entrez, SeqIO
from Bio.bgzf import BgzfWriter


def esearch_accs(db, term, retmax=1000):
    # get number of results
    with Entrez.esearch(db=db, term=term, idtype="acc", retmax=retmax, usehistory=True) as handle:
        rec = Entrez.read(handle)
    # page through results
    kwargs = dict(WebEnv=rec["WebEnv"], QueryKey=rec["QueryKey"], RetMax=rec["RetMax"] or 1)
    for i in range(0, int(rec["Count"]), retmax):
        with Entrez.esearch(db=db, term=term, idtype="acc", retstart=i, **kwargs) as handle:
            yield Entrez.read(handle)["IdList"]


def parse_argv(argv):
    parser = ArgumentParser(
        description="update a repository of indexed sequence files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("repo", type=Path, help="the target file to create/update")
    parser.add_argument("-db", default="nuccore", help="the NCBI database")
    parser.add_argument("-term", help="the NCBI query term")
    parser.add_argument("-rettype", default="fasta", help="the sequence file format")
    parser.add_argument("-retmax", type=int, default=1000, help="the records to post at a time")
    parser.add_argument("-xmdat", action="store_true", help="flag to disable modified date limit")
    parser.add_argument("-email", default="", help="the e-mail to identify yourself to NCBI")

    args = parser.parse_args(argv)

    return args


def main(argv):
    args = parse_argv(argv[1:])

    # set e-mail for identification to NCBI
    Entrez.email = args.email

    # repo directory
    makedirs(args.repo.parent, exist_ok=True)
    # repo database
    path_db = args.repo.with_suffix(".db")
    # repo log
    path_log = args.repo.with_suffix(".log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=(logging.FileHandler(path_log), logging.StreamHandler()),
    )
    logging.info(argv)

    # metadata
    accs, fdat, meta = set(), {}, dict(mdat=datetime.now().strftime("%Y/%m/%d"))
    if path_db.exists():
        with sqlite3.connect(path_db) as conn:
            # 0 -> key
            accs = {row[0] for row in conn.execute("SELECT key FROM offset_data")}
            # file_number -> name
            fdat = OrderedDict(conn.execute("SELECT * FROM file_data"))
            # key -> value
            meta = OrderedDict(conn.execute("SELECT * FROM meta_data"))
            # override args
            args.db = meta["db"]
            args.rettype = meta["format"]
            args.term = args.term if args.xmdat else f"{meta['term']} AND {meta['mdat']}:3000[MDAT]"

    # remote - local accessions
    logging.info(args.term)
    accs = list(set(chain.from_iterable(esearch_accs(args.db, args.term, args.retmax))) - accs)
    logging.info(f"count = {len(accs)}")

    paths = []
    width = len(str(args.retmax))
    for i, j in enumerate(range(0, len(accs), args.retmax), start=1):
        # fetch
        k = min(len(accs), j + args.retmax)
        csv = ",".join(accs[j : j + args.retmax])
        with Entrez.efetch(args.db, id=csv, rettype=args.rettype, retmode="text") as handle:
            path = args.repo.parent / f"{args.repo.name}-{i}.{args.rettype}.bgz.tmp"
            # compress
            with BgzfWriter(path) as stream:
                print(handle.read(), file=stream)
        paths.append(path)
        logging.info(f"{j:0{width}} - {k:0{width}} {k / len(accs) * 100:06.2f}%")

    # truthy indicates new accessions
    if paths:
        # combine previous files with new ones
        paths = list(map(Path, fdat.values())) + paths
        # rename with zero-fill
        width = len(str(len(paths)))
        paths = {
            ele: ele.with_name(f"{args.repo.name}-{idx:0{width}}.{args.rettype}.bgz")
            for idx, ele in enumerate(paths, start=1)
        }
        for key, val in paths.items():
            logging.info(f"{key} -> {val}")
            key.rename(val)
        try:
            path_tmp = path_db.with_suffix(".tmp")
            path_tmp.exists() and path_tmp.unlink()
            print("index...")
            SeqIO.index_db(str(path_tmp), list(map(str, paths.values())), args.rettype)
            # update metadata
            with sqlite3.connect(path_tmp) as conn:
                conn.execute(
                    "INSERT INTO meta_data VALUES ('db', ?), ('term', ?), ('mdat', ?)",
                    (args.db, args.term, meta["mdat"]),
                )
            path_tmp.rename(path_db)
        except Exception as e:
            logging.error(e)
            # revert original path names
            for key, val in paths.items():
                logging.info(f"{val} -> {key}")
                val.exists() and val.rename(key)

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
