#!/usr/bin/env python3

import sys
from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL


def parse_argv(argv):
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"file", type=FileType()
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.file as file:
		for record in SeqIO.parse(file, "genbank"):
			feats = getattr(record, "features", [{}])
			quals = getattr(feats[0], "qualifiers", {})
			db_xref = quals.get("db_xref", ["taxon:0"])
			taxon = next((ele for ele in db_xref if ele.startswith("taxon:")))
			taxon = taxon.split(":")[-1]
			print(record.id, taxon)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
