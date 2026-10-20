#!/usr/bin/env python3
from pathlib import Path


def gt(value):
    return "0/0" if value == 0 else "1/1"


out = Path(__file__).resolve().parent
cultivars = [f"C{i:02d}" for i in range(40)]
wild = [f"W{i:02d}" for i in range(4)]
(out / "cultivar.list").write_text("\n".join(cultivars) + "\n")
(out / "wild.list").write_text("\n".join(wild) + "\n")
with (out / "synthetic.vcf").open("w") as handle:
    handle.write("##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=2000000>\n")
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(cultivars + wild) + "\n")
    for site, pos in enumerate(range(5000, 2000000, 10000)):
        values = []
        for i in range(40):
            if pos < 1600000:
                # Stable 50/50 backbone haplotypes; a short interval has too few
                # common SNPs because records there are made rare below.
                value = 0 if i < 20 else 2
            else:
                # A sustained membership change defines another block.
                value = 0 if i % 2 == 0 else 2
            # Rare variants create recent structure without changing most common sites.
            if site % 13 == 0 and i == (site % 4):
                value = 2
            values.append(gt(value))
        wild_values = ["0/0", "0/0", "1/1", "1/1"]
        handle.write(f"Chr01\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(values + wild_values) + "\n")
