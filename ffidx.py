#!/usr/bin/env python3

import os
import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def ffidx_search(index):
	target = Path(index).expanduser().with_suffix(".idx")

	if not target.exists():
		for path in map(Path, os.environ.get("FFIDX", os.getcwd()).split(":")):
			path_index = path.expanduser().joinpath(target)
			print(path_index)
			if path_index.exists():
				target = path_index
				break

	return target


def accverize(index, keys):
	with sqlite3.connect(index) as conn:
		curs = conn.cursor()
		for key in keys:
			curs.execute(
				"""
				SELECT key
				FROM offset_data
				WHERE key LIKE ? || '.%'
				ORDER BY length(key) DESC, key DESC
				LIMIT 1
				;
				""",
				(key,)
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
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"index",
		help="the SQLite index"
	)
	parser.add_argument(
		"-filenames", "--filenames",
		nargs="+",
		help="the list of sequence files to index"
	)
	parser.add_argument(
		"-all", "-all", "-dump", "--dump",
		dest="all",
		action="store_true",
		help="the flag to dump all of the records"
	)
	parser.add_argument(
		"-descriptions", "--descriptions", "-headers", "--headers",
		action="store_true",
		help="the flag to only output the sequence descriptions"
	)
	parser.add_argument(
		"-entry", "--entry", "-accessions", "--accessions",
		nargs="+",
		help="the accessions to retrieve"
	)
	parser.add_argument(
		"-entry-batch", "--entry-batch",
		type=FileType(),
		help="the file of accessions to retrieve"
	)
	parser.add_argument(
		"-keyerror", "--keyerror",
		action="store_true",
		help="the flag to exit on key error (if the accession isn't found)"
	)
	parser.add_argument(
		"-no-version", "--no-version",
		dest="no_ver",
		action="store_true",
		help="the flag to indicate that the accessions are missing a version"
	)
	parser.add_argument(
		"-fmt-idx", "--fmt-idx",
		help="the sequence file format of the indexed files, optional if reloading"
	)
	parser.add_argument(
		"-fmt-out", "--fmt-out",
		default="fasta",
		help="the sequence file format (output)"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	index = str(ffidx_search(args.index))

	db = SeqIO.index_db(index, filenames=args.filenames, format=args.fmt_idx)

	keys = []
	if args.all:
		keys = db.keys()
	else:
		if args.entry:
			keys = args.entry
		if args.entry_batch:
			with args.entry_batch as file:
				keys += list(map(str.strip, file))

	keys = accverize(index, keys) if args.no_ver else keys
	records = keygetter(db, keys, keyerror=args.keyerror)
	if args.descriptions:
		print(*(record.description for record in records), sep="\n")
	else:
		SeqIO.write(records, sys.stdout, args.fmt_out)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
