#!/usr/bin/env python3

import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL
from subprocess import PIPE, Popen

from Bio import Entrez, SeqIO
from Bio.bgzf import BgzfWriter

ext_to_fmt = {
	"faa": "fasta",
	"fna": "fasta",
	"fas": "fasta",
	"fasta": "fasta",
	"gbk": "genbank",
	"genbank": "genbank"
}


def batchify(entries, size=10):
	batch = []
	for i, e in enumerate(entries, start=1):
		batch.append(e)
		if i % size == 0:
			yield batch
			batch = []

	if batch:
		yield batch


def parse_argv(argv):
	parser = ArgumentParser(
		description="repo auto util",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	# local
	parser.add_argument(
		"repo",
		type=Path,
		help="the target file to create/update"
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format",
		default="fasta"
	)
	parser.add_argument(
		"-redo", "--redo", action="store_true"
	)
	# remote
	parser.add_argument(
		"-db", "--db", "-database", "--database",
		default="nuccore",
		help="the database"
	)
	parser.add_argument(
		"-term",
		help="the term"
	)
	parser.add_argument(
		"-post-size", "--post-size",
		type=int,
		default=1000,
		help="the number of records to post at a time"
	)
	parser.add_argument(
		"-fetch-size", "--fetch-size",
		type=int,
		default=100,
		help="the number of records to fetch at a time"
	)
	parser.add_argument(
		"-email", "--email",
		default="",
		help="the e-mail"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	Entrez.email = args.email

	repo = args.repo
	ext = args.fmt
	fmt = ext_to_fmt[ext]

	db = args.db
	term = args.term

	path_idx = repo.with_suffix(".idx")

	meta, files, accs = {}, [], set()
	if args.redo:
		path_idx.exists() and path_idx.unlink() and print("unlink: ", str(path_idx))
		for path in path_idx.parent.glob(f"{repo}.[0-9].{ext}.bgz"):
			print("unlink: ", str(path))
			path.unlink()
	elif path_idx.exists():
		with sqlite3.connect(path_idx) as conn:
			meta = dict(conn.execute("SELECT * FROM meta_data").fetchall())
			accs = {row[0] for row in conn.execute("SELECT key FROM offset_data")}

	db = meta.get("db", db)
	fmt = meta.get("format", fmt)
	ext = meta.get("ext", ext)
	term = meta.get("term", term)

	# get the number of remote accessions
	with Entrez.esearch(db=db, term=term, idtype="acc") as file:
		record = Entrez.read(file)

	# get the remote accessions
	count = int(record["Count"])
	print("count: ", count, file=sys.stderr)

	with Entrez.esearch(db=db, term=term, idtype="acc", retmax=count) as file:
		record = Entrez.read(file)
		accs = set(record["IdList"]) - accs

	# accession diff
	print("new: ", len(accs), file=sys.stderr)

	if accs:
		# cat file | epost -db db | efetch -format fmt > target
		for batch in batchify(accs, args.post_size):
			print("download:", *batch[:5], "...", file=sys.stderr)
			kwargs = dict(stdin=PIPE, stdout=PIPE, universal_newlines=True)
			cmd1, cmd2 = ["epost", "-db", db], ["efetch", "-format", fmt]

			filenames = []
			if path_idx.exists():
				with sqlite3.connect(path_idx) as conn:
					rows = conn.execute("SELECT name FROM file_data ORDER BY file_number").fetchall()
					filenames = [repo.parent / row[0] for row in rows]

			path_rec = repo.with_suffix(f".{len(filenames)}.{ext}.bgz")
			with BgzfWriter(str(path_rec)) as file:
				with Popen(cmd1, **kwargs) as pipe1, Popen(cmd2, **kwargs) as pipe2:
					stdout, stderr = pipe1.communicate("\n".join(batch))
					stdout, stderr = pipe2.communicate(stdout)
					print(stdout, file=file)

			path_idx_tmp = path_idx.with_suffix(".tmp")
			filenames += [str(path_rec)]
			SeqIO.index_db(str(path_idx_tmp), filenames=filenames, format=fmt)
			with sqlite3.connect(path_idx_tmp) as conn:
				conn.execute(
					"INSERT INTO meta_data VALUES ('db', ?), ('term', ?), ('ext', ?)",
					(db, term, ext)
				)
			path_idx_tmp.replace(path_idx)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
