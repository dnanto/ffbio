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
	"fas": "fasta",
	"faa": "fasta",
	"fna": "fasta",
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


def get_index_paths(path_idx):
	filenames = []

	if path_idx.exists():
		with sqlite3.connect(path_idx) as conn:
			rows = conn.execute("SELECT name FROM file_data ORDER BY file_number").fetchall()
			filenames = [path_idx.parent / row[0] for row in rows]

	return filenames


def reindex(path_idx, filenames, fmt, db, term, ext):
	path_idx_tmp = path_idx.with_suffix(".tmp")
	path_idx_tmp.exists() and path_idx_tmp.unlink()
	SeqIO.index_db(str(path_idx_tmp), filenames=filenames, format=fmt)
	with sqlite3.connect(path_idx_tmp) as conn:
		conn.execute(
			"INSERT INTO meta_data VALUES ('db', ?), ('term', ?), ('ext', ?)",
			(db, term, ext)
		)
	path_idx_tmp.replace(path_idx)


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
		"-cache", "--cache", nargs="+"
	)
	parser.add_argument(
		"-cache-fmt", "--cache-fmt", default="fasta"
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

	# load cache
	cache = {}
	if args.cache:
		cache = SeqIO.index_db(":memory:", filenames=args.cache, format=ext_to_fmt[args.cache_fmt])

	# redo or load metadata
	meta, files, accs = {}, [], set()
	if args.redo:
		for path in get_index_paths(path_idx):
			print("unlink: ", str(path))
			path.unlink()
		if path_idx.exists():
			print("unlink: ", str(path_idx))
			path_idx.unlink()
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

	# accession diff length
	print("new: ", len(accs), file=sys.stderr)

	keys = accs & set(cache)
	accs -= keys
	print("cached:", len(keys), file=sys.stderr)
	print("download:", len(accs), file=sys.stderr)

	# update from cache
	if keys:
		filenames = list(map(str, get_index_paths(path_idx)))
		filenames += [str(path_idx.with_suffix(f".{len(filenames)}.{ext}.bgz"))]

		with BgzfWriter(filenames[-1]) as file:
			SeqIO.write((cache[key] for key in cache), file, fmt)

		reindex(path_idx, filenames, fmt, db, term, ext)

	if accs:
		# cat file | epost -db db | efetch -format fmt > target
		for batch in batchify(accs, args.post_size):
			print("download:", *batch[:5], "...", file=sys.stderr)
			kwargs = dict(stdin=PIPE, stdout=PIPE, universal_newlines=True)
			cmd1, cmd2 = ["epost", "-db", db], ["efetch", "-format", fmt]

			filenames = list(map(str, get_index_paths(path_idx)))
			filenames += [str(path_idx.with_suffix(f".{len(filenames)}.{ext}.bgz"))]

			with BgzfWriter(filenames[-1]) as file:
				with Popen(cmd1, **kwargs) as pipe1, Popen(cmd2, **kwargs) as pipe2:
					stdout, stderr = pipe1.communicate("\n".join(batch))
					stdout, stderr = pipe2.communicate(stdout)
					print(stdout, file=file)

			reindex(path_idx, filenames, fmt, db, term, ext)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
