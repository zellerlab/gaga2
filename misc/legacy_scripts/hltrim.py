import os
import sys
import gzip
import argparse
import pathlib
import contextlib

from collections import Counter


def read_fastq_stream(fq_in):
	for line in fq_in:
		_id, seq, _, qual = line, next(fq_in), next(fq_in), next(fq_in)
		yield _id.strip(), seq.strip(), qual.strip()


def process_fastq_record(fq_rec, cutlen):
	if fq_rec is None:
		return None

	_id, seq, qual = fq_rec
	_len = len(seq)

	if _len < cutlen:
		return None

	return _id, seq[:cutlen], qual[:cutlen]


def main():
	ap = argparse.ArgumentParser()
	# ap.add_argument("--r1", type=str)
	# ap.add_argument("--r2", type=str)
	ap.add_argument("reads", nargs="*", type=str)
	ap.add_argument("--outdir", "-o", type=str, default="output")
	ap.add_argument("--cutlen", "-c", type=str)
	args = ap.parse_args()

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)

	cutlen = args.cutlen.split(",")
	if len(cutlen) == 1:
		r1cut, r2cut = int(cutlen[0]), None
	else:
		r1cut, r2cut, *_ = map(int, cutlen)

	reads = sorted(args.reads)
	r1 = reads[0]
	r2 = reads[1] if len(reads) > 1 else None

	_open = gzip.open if r1.endswith(".gz") else open


	r1_in = _open(r1, "rt")
	r2_in = _open(r2, "rt") if r2 else contextlib.nullcontext()

	r1_out = _open(os.path.join(args.outdir, r1), "wt")
	r2_out = _open(os.path.join(args.outdir, r2), "wt") if r2 else contextlib.nullcontext()

	with r1_in, r2_in, r1_out, r2_out:
		r1_stream = read_fastq_stream(r1_in)
		r2_stream = read_fastq_stream(r2_in) if r2 else None

		while True:
			try:
				r1_rec = process_fastq_record(next(r1_stream), r1cut)
			except StopIteration:
				break
			if r2_stream:
				try:
					r2_rec = process_fastq_record(next(r2_stream), r2cut)
				except StopIteration:
					raise ValueError("r2 cannot be exhausted before r1")
			else:
				r2_rec = None

			if r1_rec is not None and (r2_rec is not None or r2_stream is None):
				print(*r1_rec[0:2], "+", r1_rec[-1], sep="\n", file=r1_out)
				if r2_rec is not None:
					print(*r2_rec[0:2], "+", r2_rec[-1], sep="\n", file=r2_out)
				
			
			
	
		

if __name__ == "__main__":
    main()
