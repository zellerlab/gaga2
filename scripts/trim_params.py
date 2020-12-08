import sys
import json
import argparse
import gzip

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("trim_params", type=str)
	ap.add_argument("reads", type=str)
	args = ap.parse_args()

	d = json.load(open(args.trim_params))

	"/g/scb/bork/schudoma/16S/testrun/input2/SRR5906274/SRR5906274_1.fastq.gz"
	fopen = gzip.open if args.reads.endswith(".gz") else open
	with fopen(args.reads, "rt") as reads_in:
		seq = [next(reads_in) for i in range(2)][1]

	print(*d[0]["trimPosition"], len(seq), sep="\t")
	

	"""
	>>> import json
	>>> d=json.load(open("/g/scb/bork/schudoma/16S/work/db/62cf129e55a61dd6b11e32d9c27e76/figaro_out/trimParameters.json"))
	>>> type(d)
	<class 'list'>
	>>> d[0]
	{'trimPosition': [285, 238], 'maxExpectedError': [4, 4], 'readRetentionPercent': 74.95, 'score': 56.95061011040093}
	>>>
	"""

if __name__ == "__main__":
	main()
