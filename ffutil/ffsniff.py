#!/usr/bin/env python3

import fileinput
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from io import StringIO, IOBase
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def openhook(fn, mode):
	return fn if isinstance(fn, IOBase) else open(fn, mode)


def ffsniff(stream, size=80):
	sample = stream.read(size).lstrip()

	fmt = None

	if sample.startswith(">"):
		fmt = "fasta"
	elif sample.startswith("LOCUS"):
		fmt = "genbank"
	else:
		# todo: other formats...
		pass

	file = fileinput.input(files=(StringIO(sample), stream), openhook=openhook)
	return file, fmt


def parse_argv(argv):
	parser = ArgumentParser(
		description="filter sequence records by header and/or length",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.file as file:
		f, fmt = ffsniff(file)
		for record in SeqIO.parse(f, fmt):
			print(record.id)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
