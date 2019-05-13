#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import OrderedDict
from itertools import groupby
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

from ffkit.ffparse import ffsniff, ffparse


def keyfunc(record):
	return seguid(record.seq)


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
		"-map", "--map", type=Path,
		help="the path prefix for the output files"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# aggregate sequences by unique hash
	with args.file as file:
		sample, fmt = ffsniff(file)
		keyfunc = lambda val: seguid(val.seq)
		records = OrderedDict(
			(key, list(val)) for key, val in
			groupby(sorted(ffparse(file, fmt, sample=sample), key=keyfunc), keyfunc)
		)

	# map unique sequence index to unique hash
	if args.map:
		with args.map.open("w") as file:
			width = len(str(len(records)))
			print("idx", "key", "id", "description", sep="\t", file=file)
			for idx, ele in enumerate(records.items(), start=1):
				idx, key, val = f"{idx:0{width}d}", ele[0], ele[1]
				for record in val:
					print(idx, key, record.id, record.description, sep="\t", file=file)
					record.id = idx

	# output unique sequences
	SeqIO.write((val[0] for val in records.values()), sys.stdout, fmt)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
