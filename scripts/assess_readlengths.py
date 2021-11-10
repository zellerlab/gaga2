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

	read_lengths = Counter()
	for fastq_report in glob.glob(os.path.join(args.input_dir, "*fastqc_data.txt")):
		read_lengths.update(parse_fastqc_report(fastq_report))

	yields = list()
	for length, count in read_lengths.items():
		yields.append((length, sum(v for k, v in read_lengths.items() if k >= length), sum(v * length for k, v in read_lengths.items() if k >= length)))

	for length, reads, bases in sorted(yields, key=lambda x:(x[2], x[1]), reverse=True):
		print(length, reads, bases, sep="\t")
	
		

	# print(*read_lengths.items(), sep="\n")
		
	
	# print(read_lengths)
	


if __name__ == "__main__":
	main()
