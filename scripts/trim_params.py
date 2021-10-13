import sys
import json
import argparse
import gzip

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("trim_params", type=str)
	args = ap.parse_args()

	# {'trimPosition': [285, 238], 'maxExpectedError': [4, 4], 'readRetentionPercent': 74.95, 'score': 56.95061011040093}
	try:
		d = json.load(open(args.trim_params))
		print(*d[0]["trimPosition"])
	except:
		print("-1 -1")


if __name__ == "__main__":
	main()
