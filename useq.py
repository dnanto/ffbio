#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import OrderedDict
from itertools import groupby
from os.path import splitext
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute the unique set of sequences",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"-fmt", "--fmt", default="fasta"
	)
	parser.add_argument(
		"-echo", "--echo", action="store_true"
	)
	parser.add_argument(
		"-pfx", "--pfx", "-prefix", "--prefix",
		help="the path prefix for the output files"
	)

	args = parser.parse_args(argv)
	args.pfx = args.pfx or args.file.name
	args.echo = args.file.name == "<stdin>"

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# aggregate sequences by unique hash
	with args.file as file:
		records = OrderedDict(
			(key, list(val)) for key, val in
			groupby(SeqIO.parse(file, args.fmt), key=lambda val: seguid(val.seq))
		)

	# map unique sequence index to unique hash
	idx_to_key = OrderedDict()
	with open(args.pfx + ".tsv", "w") as file:
		print("idx", "key", "id", "description", sep="\t", file=file)
		idx_width = len(str(len(records)))
		for idx, ele in enumerate(records.items(), start=1):
			idx, key, val = f"{idx:{idx_width}d}", ele[0], ele[1]
			idx_to_key[idx] = key
			for record in val:
				print(idx, key, record.id, record.description, sep="\t", file=file)
				record.id, record.description = idx, ""

	# output unique sequences
	fname = args.file.name
	file = sys.stdout if args.echo else (args.pfx + "." + splitext(fname)[-1].lstrip("."))
	SeqIO.write((records[val][0] for val in idx_to_key.values()), file, args.fmt)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
