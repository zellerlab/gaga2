import os
import sys
import gzip
import argparse
import pathlib

from collections import Counter

FASTQ_ENDINGS = (".fq", ".fastq", ".fq.gz", ".fastq.gz")

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

    _, _, files = next(os.walk(args.indir))
    files = iter(sorted(os.path.join(args.indir, f) for f in files if any(f.endswith(suffix) for suffix in FASTQ_ENDINGS)))
    for r1 in files:
        r2 = next(files)
        r1_qdips, len_r1 = analyse_file(r1)
        r2_qdips, len_r2 = analyse_file(r2)
        qdips.update(r1_qdips)
        qdips.update(r2_qdips)
        read_files.append((r1, r2))

    figaro_path = os.path.join(args.outdir, "figaro")
    dada2_path = os.path.join(args.outdir, "dada2")
    pathlib.Path(figaro_path).mkdir(exist_ok=True, parents=True)
    pathlib.Path(dada2_path).mkdir(exist_ok=True, parents=True)
    print(*sorted(qdips.items()), sep="\n")
    cutpos = max(c for c, count in qdips.items() if c <= args.max_pos and (count > 1 or c == -1))
    print(read_files)
    print(f"ltrim: cut at {cutpos}, lengths={lengths}")
    for r1, r2 in read_files:
        # path = os.path.join(dada2_path, os.path.basename(os.path.dirname(r1)))
        path = os.path.join(args.outdir, os.path.basename(os.path.dirname(r1)))
        pathlib.Path(path).mkdir(exist_ok=True, parents=True)
        trim_file(r1, path, new_start=cutpos + 1)
        trim_file(r2, path, new_start=cutpos + 1)
        # r1_base, r2_base = map(os.path.basename, (r1, r2))
        #  try:
        #    os.symlink(os.path.join(path, r1_base), os.path.join(figaro_path, r1_base))
        #except:
        #     pass
        #try:
        #    os.symlink(os.path.join(path, r2_base), os.path.join(figaro_path, r2_base))
        #except:
        #    pass


if __name__ == "__main__":
    main()
