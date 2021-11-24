import argparse
import glob
import os

from collections import Counter

def parse_fastqc_report(f):
	readlengths = list()
	active = False
	for line in open(f):
		active = active or (not active and line.startswith(">>Sequence Length Distribution"))
		if active:
			if line.startswith(">>END_MODULE"):
				break
			if line[0] not in (">", "#"):
				length, count = line.split("\t")
				readlengths.append((int(length.split("-")[0]), float(count)))
	return dict(item for item in readlengths if item[1] > 0)


def parse_bbduk_hist(f):
	return {
		int(line.strip().split("\t")[0]): int(line.strip().split("\t")[1])
		for line in open(f)
		if line[0] != "#"
	}


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir", type=str, default=".")
	ap.add_argument("--amplicon_length", type=int, required=True)
	ap.add_argument("--min_overlap", type=int, required=True)
	args = ap.parse_args()

	minlength = 0.5 * args.amplicon_length + args.min_overlap
	# minlength = 0

	read_yields = {}
	for mate in (1, 2):
		read_lengths = Counter()
		#for fastq_report in glob.glob(os.path.join(args.input_dir, f"*{mate}_fastqc_data.txt")):
		#	read_lengths.update(parse_fastqc_report(fastq_report))
		for hist in glob.glob(os.path.join(args.input_dir, f"*R{mate}.post_lhist.txt")):
			read_lengths.update(parse_bbduk_hist(hist))

		yields = list()
		for length, count in read_lengths.items():
			yields.append(
				(length, sum(v for k, v in read_lengths.items() if k >= length), sum(v * length for k, v in read_lengths.items() if k >= length))
			)
		for length, reads, bases in sorted(yields, key=lambda x:(x[2], x[1]), reverse=True):
			read_yields.setdefault(mate, list()).append((length, reads, bases))

	r1_lengths = read_yields.get(1, list())
	r2_lengths = read_yields.get(2, list())

	if not r2_lengths:
		all_lengths = r1_lengths
		is_hom = len(r1_lengths) == 1
	else:
		is_hom = len(r1_lengths) == len(r2_lengths) == 1
		smaller = min(len(r1_lengths), len(r2_lengths))
		all_lengths = [x + y for x, y in zip(r1_lengths[:smaller], r2_lengths[:smaller])]

	for item in all_lengths:
		if len(item) == 3 or item[0] + item[3] > minlength:
			print(*item, sep="\t")

	if is_hom:
		open("READSET_HOMOGENEOUS", "wt").close()
		


if __name__ == "__main__":
	main()
