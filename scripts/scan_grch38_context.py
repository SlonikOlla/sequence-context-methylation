#!/usr/bin/env python3
import gzip, csv, os, time
from collections import Counter, defaultdict

BASES = b"ACGT"
VALID = set(BASES)

def clean_path(s):
    return s.strip().strip('"').strip("'")

def open_fasta(path):
    return gzip.open(path, "rb") if path.lower().endswith(".gz") else open(path, "rb")

def write_counter(path, counter, headers=("motif","count")):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for k, v in sorted(counter.items()):
            w.writerow([k, v])

def write_pos(path, pos_counts):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group","position","A","C","G","T","total"])
        for group in sorted(pos_counts):
            for pos in sorted(pos_counts[group]):
                c = pos_counts[group][pos]
                w.writerow([group, pos, c["A"], c["C"], c["G"], c["T"], sum(c[b] for b in "ACGT")])

def add_pos(pos_counts, group, pos, base):
    if base in "ACGT":
        pos_counts[group][pos][base] += 1

def process_contig(name, seq, cpg6, cpg8, cph7, pos_cpg, pos_cph, chrom_rows):
    n = len(seq)
    counts = Counter()
    i = seq.find(b"C")
    while i != -1:
        if i + 1 < n:
            nxt = chr(seq[i+1])
            if nxt == "G":
                counts["CpG"] += 1
                # NNCGNN: two bases before C, CG, two bases after G.
                if i >= 2 and i + 3 < n:
                    s6 = seq[i-2:i+4]
                    if all(x in VALID for x in s6):
                        m6 = s6.decode("ascii")
                        cpg6[m6] += 1
                        if m6 == "TACGTA": counts["TACGTA"] += 1
                        if m6 == "TCCGTA": counts["TCCGTA"] += 1
                        if m6 == "TACGGA": counts["TACGGA"] += 1
                # NNNCGNNN: three bases before C, CG, three bases after G.
                if i >= 3 and i + 4 < n:
                    s8 = seq[i-3:i+5]
                    if all(x in VALID for x in s8):
                        cpg8[s8.decode("ascii")] += 1
                        # For CpG, +1..+3 are defined as bases after the G,
                        # so the invariant G itself is not treated as a predictor.
                        for p, j in [(-3,i-3),(-2,i-2),(-1,i-1),(1,i+2),(2,i+3),(3,i+4)]:
                            add_pos(pos_cpg, "CpG", p, chr(seq[j]))
            elif nxt in "ACT":
                group = "Cp" + nxt
                counts[group] += 1
                # CpH 7-mer: NNNCNNN, central C with 3 bases on each side.
                if i >= 3 and i + 3 < n:
                    s7 = seq[i-3:i+4]
                    if all(x in VALID for x in s7):
                        cph7[group][s7.decode("ascii")] += 1
                        for p, j in [(-3,i-3),(-2,i-2),(-1,i-1),(1,i+1),(2,i+2),(3,i+3)]:
                            add_pos(pos_cph, group, p, chr(seq[j]))
        i = seq.find(b"C", i+1)

    chrom_rows.append([name, n, counts["CpG"], counts["CpA"], counts["CpC"], counts["CpT"],
                       counts["TACGTA"], counts["TCCGTA"], counts["TACGGA"]])

def main():
    print("GRCh38 sequence-context scanner")
    print("--------------------------------")
    path = clean_path(input("Paste the full path to your GENCODE GRCh38 primary-assembly FASTA (.fa.gz):\n> "))
    if not os.path.isfile(path):
        print("\nERROR: File not found:", path)
        input("Press Enter to close...")
        return 1

    outdir = os.path.join(os.path.dirname(path), "GRCh38_context_results")
    os.makedirs(outdir, exist_ok=True)

    cpg6 = Counter()
    cpg8 = Counter()
    cph7 = {"CpA": Counter(), "CpC": Counter(), "CpT": Counter()}
    pos_cpg = defaultdict(lambda: defaultdict(Counter))
    pos_cph = defaultdict(lambda: defaultdict(Counter))
    chrom_rows = []

    start = time.time()
    name = None
    seq = bytearray()
    contigs = 0

    def flush():
        nonlocal name, seq, contigs
        if name is None:
            return
        contigs += 1
        print(f"Scanning {name} ({len(seq):,} bp)...", flush=True)
        process_contig(name, bytes(seq), cpg6, cpg8, cph7, pos_cpg, pos_cph, chrom_rows)
        seq = bytearray()

    with open_fasta(path) as f:
        for raw in f:
            if raw.startswith(b">"):
                flush()
                name = raw[1:].split()[0].decode("utf-8", errors="replace")
            else:
                seq.extend(raw.strip().upper())
        flush()

    write_counter(os.path.join(outdir, "CpG_NNCGNN_counts.csv"), cpg6)
    write_counter(os.path.join(outdir, "CpG_NNNCGNNN_counts.csv"), cpg8)
    for group in ("CpA","CpC","CpT"):
        write_counter(os.path.join(outdir, f"{group}_NNNCNNN_counts.csv"), cph7[group])
    write_pos(os.path.join(outdir, "CpG_flanking_position_counts.csv"), pos_cpg)
    write_pos(os.path.join(outdir, "CpH_flanking_position_counts.csv"), pos_cph)

    with open(os.path.join(outdir, "chromosome_summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["contig","length_bp","CpG","CpA","CpC","CpT","TACGTA","TCCGTA","TACGGA"])
        w.writerows(chrom_rows)

    with open(os.path.join(outdir, "README.txt"), "w", encoding="utf-8") as f:
        f.write(
            "GRCh38 sequence-context scan\n"
            "CpG_NNCGNN_counts.csv: 2 bases before C, CG, 2 bases after G.\n"
            "CpG_NNNCGNNN_counts.csv: 3 bases before C, CG, 3 bases after G.\n"
            "CpG_flanking_position_counts.csv: positions -3,-2,-1 before C and +1,+2,+3 after G.\n"
            "CpA/CpC/CpT_NNNCNNN_counts.csv: central C with 3 bases on each side; +1 is A/C/T.\n"
            "Only A/C/G/T-only contexts are counted in motif/position tables; ambiguous bases are excluded.\n"
            "Primary assembly FASTA may include non-chromosomal primary scaffolds; chromosome_summary.csv preserves contig identity.\n"
        )

    elapsed = (time.time() - start) / 60
    print("\nDONE.")
    print("Results folder:", outdir)
    print(f"Processed {contigs} contigs in {elapsed:.1f} minutes.")
    print("\nPlease ZIP the GRCh38_context_results folder and upload that ZIP to ChatGPT.")
    input("\nPress Enter to close...")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
