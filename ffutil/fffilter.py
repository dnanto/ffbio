#!/usr/bin/env python3

import operator as op
import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO

ops = {
	"<": op.lt,
	"<=": op.le,
	"=": op.eq,
	"==": op.eq,
	"!=": op.eq,
	">=": op.ge,
	">": op.gt,
	"in": op.contains,
}


def parse_query(query):
	cmd = None

	tokens = query.strip().split(" ", maxsplit=1)

	if re.match(r"\s*N?LEN", tokens[0], re.I):
		match = re.search(r"([><=]+|in)\s*((\d+)\s*-\s*(\d+)|\d+)", tokens[1], re.I)
		opkey = match.group(1).lower()
		if opkey == "in":
			rng = range(int(match.group(3)), int(match.group(4)) + 1)
			cmd = lambda rec: ops[opkey](rng, len(rec))
		else:
			cmd = lambda rec: ops[opkey](len(rec), int(match.group(2)))

	if re.match(r"\s*N?LIKE", tokens[0], re.I):
		tokens = tokens[1].strip().split(" ", maxsplit=1)
		cmd = lambda rec: re.search(tokens[1].strip(), getattr(rec, tokens[0].strip()))

	return (lambda x: not cmd(x)) if query.strip()[0] in "Nn" else cmd


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
	parser.add_argument(
		"query",
		help="the query to filter records"
	)
	parser.add_argument(
		"-fmt-i", "--fmt-i",
		default="fasta",
		help="the sequence file format (input)"
	)
	parser.add_argument(
		"-fmt-o", "--fmt-o",
		default="fasta",
		help="the sequence file format (output)"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	cmd = parse_query(args.query)

	with args.file as file:
		SeqIO.write(filter(cmd, SeqIO.parse(file, args.fmt_i)), sys.stdout, args.fmt_o)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
