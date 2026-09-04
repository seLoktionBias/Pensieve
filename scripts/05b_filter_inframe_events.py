#!/usr/bin/env python3
"""Write a copy of 03_<GENE>.alignment_events.tsv without in-frame indels.

Drop rule is exactly 05_plot_events.R's own is_inframe_indel test:
    character_class == "indel"  AND  as.numeric(length_mod_3) == 0

That covers both categories of non-disrupting indel: shared ancestral ones
(marker classes shared_inframe_insertion / shared_inframe_deletion) and
lineage-specific ones (marker class "other"). Frameshift indels, terminal
indels and every stop_mask event are passed through byte-for-byte, so the
second figure carries only frameshift indels and pseudogenizing events.
"""
import sys

src, dst = sys.argv[1], sys.argv[2]
raw = open(src, "rb").read().decode()
eol = "\r\n" if "\r\n" in raw else "\n"
lines = raw.splitlines()
header = lines[0].split("\t")
i_class, i_mod3 = header.index("character_class"), header.index("length_mod_3")

def is_inframe(line):
    f = line.split("\t")
    if f[i_class] != "indel":
        return False
    try:
        return float(f[i_mod3]) == 0
    except ValueError:
        return False

kept = [l for l in lines[1:] if not is_inframe(l)]
with open(dst, "w", newline="") as fh:
    fh.write(lines[0] + eol)
    for l in kept:
        fh.write(l + eol)
n = len(lines) - 1
print(f"[plot] in-frame indels suppressed: {n} events -> {len(kept)} kept ({n-len(kept)} dropped)")
