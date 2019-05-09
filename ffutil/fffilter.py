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


def getele(rec, epath):
	obj = rec

	for token in epath.split("."):
		if token.isdigit():
			obj = obj[int(token)]
		elif hasattr(obj, token):
			obj = getattr(obj, token)
		elif token in obj:
			obj = obj[token]

	return obj


def parse_cmd(query):
	cmd = None

	tokens = query.strip().split(" ", maxsplit=1)

	keyword = tokens[0].strip()

	if re.match(r"^N?LEN$", keyword, re.I):
		match = re.search(r"([><=]+|in)\s*((\d+)\s*-\s*(\d+)|\d+)", tokens[1], re.I)
		opkey = match.group(1).lower()
		if opkey == "in":
			ran = range(int(match.group(3)), int(match.group(4)) + 1)
			cmd = lambda rec: ops[opkey](ran, len(rec))
		else:
			cmd = lambda rec: ops[opkey](len(rec), int(match.group(2)))

	if re.match(r"^N?LIKE(IN)?$", keyword, re.I):
		tokens = tokens[1].strip().split(" ", maxsplit=1)
		if keyword.endswith("IN"):
			cmd = lambda rec: any(
				re.search(tokens[1].strip(), str(ele), re.I)
				for ele in getele(rec, str(tokens[0].strip()))
			)
		else:
			cmd = lambda rec: (
				re.search(
					tokens[1].strip(),
					getele(rec, str(tokens[0].strip())),
					re.I
				)
			)

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
		"-fmt-in", "--fmt-in",
		default="fasta",
		help="the sequence file format (input)"
	)
	parser.add_argument(
		"-fmt-out", "--fmt-out",
		help="the sequence file format (output)"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	cmd = parse_cmd(args.query)

	fmt_in = args.fmt_in
	fmt_out = args.fmt_out or args.fmt_in

	with args.file as file:
		SeqIO.write(filter(cmd, SeqIO.parse(file, fmt_in)), sys.stdout, fmt_out)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
