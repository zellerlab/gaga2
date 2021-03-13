import os
import sys
import gzip
import argparse
import pathlib

from collections import Counter


def find_qdips(quals, phred=33):
    def is_qdip(c, phred=33):
        return ord(c) - phred <= 2
    return [i for i, q in enumerate(quals) if is_qdip(q, phred=phred)]

def analyse_file(fn, phred=33):
    f_open = gzip.open if fn.endswith("gz") else open
    qdips = Counter()
    for i, line in enumerate(f_open(fn, "rt"), start=1):
        if i % 4 == 0:
            qdips.update(find_qdips(line.strip()))
    return qdips, len(line.strip())

def trim_file(fn, outdir, new_start=0):
    target = os.path.join(outdir, os.path.basename(fn))
    if new_start:
        f_open = gzip.open if fn.endswith("gz") else open
        with f_open(fn, "rt") as fastq_in, f_open(target, "wt") as fastq_out:
            for i, line in enumerate(fastq_in, start=1):
                if i % 2 == 0:
                    line = line[new_start:]
                print(line, end="", file=fastq_out)
    else:
        try:
            os.symlink(os.path.abspath(fn), target)
        except:
            raise ValueError(f"Could not create symlink. {fn} -> {target}.")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str)
    ap.add_argument("outdir", type=str)
    ap.add_argument("--max_pos", type=int, default=10)
    args = ap.parse_args()

    qdips = Counter([-1])
    lengths = set()
    read_files = list()
    for cwd, dirs, files in os.walk(args.indir):
        if len(files) == 2:
            r1, r2 = sorted(os.path.join(cwd, f) for f in files)
            r1_qdips, len_r1 = analyse_file(r1)
            r2_qdips, len_r2 = analyse_file(r2)
            lengths.update((len_r1, len_r2))
            qdips.update(r1_qdips)
            qdips.update(r2_qdips)
            read_files.append((r1, r2))

        
    pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
    cutpos = max(c for c in qdips if c <= args.max_pos)
    # print(qdips, cutpos, lengths)
    print(f"ltrim: cut at {cutpos}, lengths={lengths}")
    for r1, r2 in read_files:
        path = os.path.join(args.outdir, os.path.basename(os.path.dirname(r1)))
        pathlib.Path(path).mkdir(exist_ok=True, parents=True)
        trim_file(r1, path, new_start=cutpos + 1)
        trim_file(r2, path, new_start=cutpos + 1)
     


if __name__ == "__main__":
    main()
