import sys
import os
import glob
import argparse
import pathlib

FASTQ_ENDINGS = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir")
	ap.add_argument("link_dir")
	args = ap.parse_args()

	args.input_dir = os.path.abspath(args.input_dir)
	pathlib.Path(args.link_dir).mkdir(parents=True, exist_ok=True)

	print("INPUTDIR", args.input_dir)
	print(os.listdir(args.input_dir))

	for wd, dirs, files in os.walk(args.input_dir, followlinks=True):
		print(wd, dirs, files)
		for f in files:
			if any(map(lambda x:f.endswith(x), FASTQ_ENDINGS)):
				os.symlink(os.path.join(wd, f), os.path.join(args.link_dir, os.path.basename(f)))
				



if __name__ == "__main__":
	main()
