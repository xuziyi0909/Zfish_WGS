"""Split a genome (.fai) into N BED files of roughly equal size.

Large chromosomes are cut into pieces; small scaffolds are packed together.
"""

fai = snakemake.input[0]
n = snakemake.params.n
outdir = snakemake.params.outdir

contigs = []
with open(fai) as fh:
    for line in fh:
        name, length = line.split("\t")[:2]
        contigs.append((name, int(length)))

total = sum(length for _, length in contigs)
target = total / n

chunks = [[] for _ in range(n)]
idx, filled = 0, 0
for name, length in contigs:
    start = 0
    while start < length:
        room = max(int(target - filled), 1)
        end = min(length, start + room) if idx < n - 1 else length
        chunks[idx].append((name, start, end))
        filled += end - start
        start = end
        if filled >= target and idx < n - 1:
            idx, filled = idx + 1, 0

for i, regions in enumerate(chunks):
    with open(f"{outdir}/chunk_{i:02d}.bed", "w") as out:
        for name, start, end in regions:
            out.write(f"{name}\t{start}\t{end}\n")
