import sys
import json
import argparse
import gzip

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("trim_params", type=str)
	#ap.add_argument("reads", type=str)
	args = ap.parse_args()

	# {'trimPosition': [285, 238], 'maxExpectedError': [4, 4], 'readRetentionPercent': 74.95, 'score': 56.95061011040093}
	try:
		d = json.load(open(args.trim_params))
		print(*d[0]["trimPosition"])
	except:
		print("-1 -1")

	#try:
	#	reads_in = gzip.open(args.reads, "rt")
	#except:
	#	reads_in = open(args.reads, "rt")
	#fopen = gzip.open if args.reads.endswith(".gz") else open
	#with fopen(args.reads, "rt") as reads_in:
	#with reads_in:
	#	seq = [next(reads_in) for i in range(2)][1]
	# print(*d[0]["trimPosition"], len(seq), sep="\t")


if __name__ == "__main__":
	main()
