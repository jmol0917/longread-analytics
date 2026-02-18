#!/usr/bin/env python3

import re
import pandas as pd
import matplotlib.pyplot as plt
import os

VCF = "results/sv/HG002_subset.sv.vcf"
OUT_TABLE = "results/tables/HG002_subset.sv.table.csv"
OUT_TSV = "results/tables/HG002_subset.sv.table.tsv"
OUT_FIG = "results/figures/HG002_subset.svlen_hist.png"
OUT_DAY4_TSV = "results/sv/summary.tsv"
OUT_DAY4_FIG = "results/sv/sv_summary.png"

os.makedirs("results/sv", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

rows = []

with open(VCF) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
        pos = int(pos)

        info_dict = {}
        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v
            else:
                info_dict[item] = True

        svtype = info_dict.get("SVTYPE", "")
        end = int(info_dict.get("END", pos))
        svlen = int(info_dict.get("SVLEN", "0"))
        support = int(info_dict.get("SUPPORT", "0"))
        vaf = float(info_dict.get("VAF", "nan"))

        rows.append({
            "chrom": chrom,
            "pos": pos,
            "end": end,
            "svtype": svtype,
            "svlen": svlen,
            "abs_svlen": abs(svlen),
            "filter": flt,
            "support": support,
            "vaf": vaf,
            "id": vid
        })

df = pd.DataFrame(rows)

df.to_csv(OUT_TABLE, index=False)
df.to_csv(OUT_TSV, sep="\t", index=False)

# Summary stats
df_pass = df[df["filter"] == "PASS"]

expected_types = ["DEL", "INS", "DUP", "INV", "BND"]

total = len(df)
pass_n = int((df["filter"] == "PASS").sum())
nonpass_n = total - pass_n

type_counts = {t: int((df["svtype"] == t).sum()) for t in expected_types}

absx = df["abs_svlen"].astype(float)
min_len = float(absx.min()) if total else float("nan")
median_len = float(absx.median()) if total else float("nan")
mean_len = float(absx.mean()) if total else float("nan")
max_len = float(absx.max()) if total else float("nan")

with open(OUT_DAY4_TSV, "w") as out:
    out.write("section\tmetric\tvalue\n")
    out.write(f"calls\ttotal\t{total}\n")
    out.write(f"calls\tPASS\t{pass_n}\n")
    out.write(f"calls\tnonPASS\t{nonpass_n}\n")
    for t in expected_types:
        out.write(f"svtype\t{t}\t{type_counts[t]}\n")
    out.write(f"size_abs_bp\tmin\t{min_len}\n")
    out.write(f"size_abs_bp\tmedian\t{median_len}\n")
    out.write(f"size_abs_bp\tmean\t{mean_len}\n")
    out.write(f"size_abs_bp\tmax\t{max_len}\n")

print("\n=== SV SUMMARY ===")
print("Total SVs:", len(df))
print("PASS SVs:", len(df_pass))
print("\nSVTYPE counts:")
print(df["svtype"].value_counts())

print("\nSize statistics (absolute):")
print(df["abs_svlen"].describe())

# Plot
plt.figure(figsize=(7,5))
plt.hist(df_pass["abs_svlen"], bins=30)
plt.xlabel("Absolute SV length (bp)")
plt.ylabel("Count")
plt.title("HG002 subset: SV size distribution (PASS)")
plt.tight_layout()
plt.savefig(OUT_FIG, dpi=200)
plt.savefig(OUT_DAY4_FIG, dpi=200)

print("\nWrote:")
print(" -", OUT_TABLE)
print(" -", OUT_TSV)
print(" -", OUT_FIG)
print(" -", OUT_DAY4_TSV)
print(" -", OUT_DAY4_FIG)

