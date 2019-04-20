#!/usr/bin/env python3

import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


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
			yield curs.fetchone()[0]
		curs.close()


def parse_argv(argv):
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"-index",
		default=":memory:",
		help="here to store the SQLite index"
	)
	parser.add_argument(
		"-filenames", "--filenames",
		nargs="+",
		help="list of strings specifying file(s) to be indexed"
	)
	parser.add_argument(
		"-accessions", "--accessions", dest="keys", nargs="+"
	)
	parser.add_argument(
		"-no-version", "--no-version", action="store_true"
	)
	parser.add_argument(
		"-fmt-in", "--fmt-in", default="fasta"
	)
	parser.add_argument(
		"-fmt-out", "--fmt-out", default="fasta"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	records = SeqIO.index_db(args.index, filenames=args.filenames, format=args.fmt_in)

	args.keys = list(accverize(args.index, args.keys)) if args.keys and args.no_version else args.keys

	if args.keys:
		SeqIO.write((records[key] for key in args.keys), sys.stdout, args.fmt_out)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
