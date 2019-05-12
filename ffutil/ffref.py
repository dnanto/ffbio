#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute the length of each record in the file",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file", type=FileType(),
		help="the sequence file"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.file as file:
		print(
			"id", "description", "length", "date",
			"title", "journal", "medline_id", "pubmed_id", "comment",
			"collection_date", "note",
			sep="\t"
		)
		for rec in SeqIO.parse(file, "genbank"):
			qual = rec.features[0].qualifiers
			collection_date = qual.get("collection_date", [""])[0]
			note = qual.get("note", [""])[0]
			for ref in rec.annotations.get("references", []):
				print(
					rec.id, rec.description, len(rec), rec.annotations["date"],
					ref.title, ref.journal, ref.medline_id, ref.pubmed_id, ref.comment,
					collection_date, note,
					sep="\t"
				)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
