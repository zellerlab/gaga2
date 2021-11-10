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


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir", type=str, default=".")
	args = ap.parse_args()

	read_yields = {}
	for r in (1, 2):
		read_lengths = Counter()
		for fastq_report in glob.glob(os.path.join(args.input_dir, f"*{r}_fastqc_data.txt")):
			read_lengths.update(parse_fastqc_report(fastq_report))

		yields = list()
		for length, count in read_lengths.items():
			yields.append(
				(length, sum(v for k, v in read_lengths.items() if k >= length), sum(v * length for k, v in read_lengths.items() if k >= length))
			)
		for length, reads, bases in sorted(yields, key=lambda x:(x[2], x[1]), reverse=True):
			read_yields.setdefault(r, list()).append((length, reads, bases))

	r1_lengths = read_yields.get(1, list())
	r2_lengths = read_yields.get(2, list())

	if not r2_lengths:
		all_lengths = r1_lengths
		is_hom = len(r1_lengths) == 1
	else:
		if len(r1_lengths) < len(r2_lengths):
			r1_lengths.extend(("NA", "NA", "NA") for i in range(len(r2_lengths) - len(r1_lengths)))
		elif len(r2_lengths) < len(r1_lengths):
			r2_lengths.extend(("NA", "NA", "NA") for i in range(len(r1_lengths) - len(r2_lengths)))
		all_lengths = [x + y for x, y in zip(r1_lengths, r2_lengths)]
		is_hom = len(r1_lengths) == len(r2_lengths) == 1

	for item in all_lengths:
		print(*item, sep="\t")

	if is_hom:
		open("READSET_HOMOGENEOUS", "wt").close()
		


if __name__ == "__main__":
	main()
