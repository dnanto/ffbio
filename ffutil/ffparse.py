#!/usr/bin/env python3

import fileinput
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from io import StringIO, IOBase
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def openhook(fn, mode):
	return fn if isinstance(fn, IOBase) else open(fn, mode)


def ffparse(handle, fmt=None, size=None):
	if fmt:
		yield from SeqIO.parse(handle, fmt)
	else:
		sample = handle.read(size).lstrip()

		fmt = None

		if sample.startswith(">"):
			fmt = "fasta"
		elif sample.startswith("LOCUS"):
			fmt = "genbank"
		else:
			# todo: other formats...
			pass

		with fileinput.input(files=(StringIO(sample), handle), openhook=openhook) as file:
			yield from SeqIO.parse(file, fmt)


def parse_argv(argv):
	parser = ArgumentParser(
		description="parse sequence records, sniff format automatically",
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
		for record in ffparse(file):
			print(record.id)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
