#!/usr/bin/env python3

import fileinput
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from io import IOBase, StringIO
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def ffsample(handle, n=1):
	try:
		for _ in range(n):
			yield next(handle)
	except StopIteration:
		pass


def ffsniff(handle, n=1):
	sample = "".join(ffsample(handle, n))

	fmt = None

	if sample.startswith(">"):
		fmt = "fasta"
	elif sample.startswith("LOCUS"):
		fmt = "genbank"
	else:
		# todo: other formats...
		pass

	return sample, fmt


def ffparse(handle, format, alphabet=None, sample=None):
	sample = StringIO(sample)
	with fileinput.input((sample, handle), openhook=openhook) as file:
		yield from SeqIO.parse(file, format, alphabet)


def openhook(fn, mode):
	return fn if isinstance(fn, IOBase) else open(fn, mode)


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
		sample, fmt = ffsniff(file)
		for record in ffparse(file, fmt, sample=sample):
			print(record.id)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
