import sys
import os
import gzip
import argparse
import pathlib

FASTQ_ENDINGS = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def check_readlengths(fn):
	lengths = set()
	fopen = gzip.open if fn.endswith(".gz") else open
	with fopen(fn, "rt") as reads_in:
		for i, line in enumerate(reads_in):
			if i % 4 == 1:
				lengths.add(len(line.strip()))
				if len(lengths) > 1:
					return False
	return True


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir", type=str)
	ap.add_argument("output_dir", type=str)
	args = ap.parse_args()


	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

	for wd, dirs, files in os.walk(args.input_dir):
		for fn in files:
			fn = os.path.join(wd, fn)
			if any(map(lambda x:fn.endswith(x), FASTQ_ENDINGS)) and not check_readlengths(fn):
				print("WARNING: read lengths are not homogenous in {fq}. Figaro will not be executed.".format(fq=fn))
				open(os.path.join(args.output_dir, "SKIP_FIGARO"), "wt").close()
				return None
	open(os.path.join(args.output_dir, "RUN_FIGARO"), "wt").close()


	


if __name__ == "__main__":
	main()	
