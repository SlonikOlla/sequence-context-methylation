#!/usr/bin/env python3
import gzip, csv, os, sys, time
from pathlib import Path
from collections import defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"

def clean_path(s):
    return s.strip().strip('"').strip("'")

def open_text(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="ignore")
    return open(path, "rt", encoding="ascii", errors="ignore")

def revcomp(s):
    return s.translate(COMP)[::-1]

def load_reference(path):
    print("\nLoading GRCh38 reference into memory (this can take several minutes)...")
    seqs, name, chunks = {}, None, []
    with open_text(path) as f:
        for line in f:
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                    print(f"  {name}: {len(seqs[name]):,} bp")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks).upper()
            print(f"  {name}: {len(seqs[name]):,} bp")
    print(f"Loaded {len(seqs)} reference sequences.")
    return seqs

def band_for(cov):
    if cov < 10:
        return None
    if cov <= 14:
        return "10-14x"
    if cov <= 19:
        return "15-19x"
    if cov <= 29:
        return "20-29x"
    return "30plus"

# Aggregate: total coverage, methylated-read equivalents, number of sites,
# and sum of site methylation percentages. This lets us report both
# coverage-weighted and site-unweighted means.
def add(agg, key, cov, pct):
    a = agg[key]
    a[0] += cov
    a[1] += cov * pct / 100.0
    a[2] += 1
    a[3] += pct

def write_agg(path, header, agg, key_unpack):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header + [
            "site_count", "total_coverage",
            "mean_methylation_percent",
            "coverage_weighted_methylation_percent"
        ])
        for key in sorted(agg, key=lambda x: str(x)):
            total_cov, meth_reads, n, sum_pct = agg[key]
            site_mean = sum_pct / n if n else float("nan")
            weighted = 100.0 * meth_reads / total_cov if total_cov else float("nan")
            vals = key_unpack(key)
            w.writerow(list(vals) + [
                n, int(total_cov), f"{site_mean:.8f}", f"{weighted:.8f}"
            ])

def process_bed(bed_path, label, context, seqs, outdir):
    context = context.upper()
    print(f"\n=== {label}: {context} ===")
    print(f"Input: {bed_path}")
    print("Minimum coverage: 10x")
    print("Large CHH files can take a while. Progress prints every 5,000,000 rows.")

    overall = defaultdict(lambda:[0.0,0.0,0,0.0])
    subtype = defaultdict(lambda:[0.0,0.0,0,0.0])
    position = defaultdict(lambda:[0.0,0.0,0,0.0])
    motif = defaultdict(lambda:[0.0,0.0,0,0.0])
    bin_subtype = defaultdict(lambda:[0.0,0.0,0,0.0])
    bin_position = defaultdict(lambda:[0.0,0.0,0,0.0])
    bin_motif = defaultdict(lambda:[0.0,0.0,0,0.0])
    qc = defaultdict(int)

    t0 = time.time()
    with open_text(bed_path) as f:
        for i, line in enumerate(f, 1):
            if i % 5_000_000 == 0:
                dt = (time.time() - t0) / 60
                print(f"  {i:,} rows scanned; {qc['kept']:,} sites kept; {dt:.1f} min")
            if not line or line[0] in "#tb":  # comments / track / browser
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                qc["short_rows"] += 1
                continue

            try:
                chrom = parts[0]
                start = int(parts[1])       # BED: 0-based
                strand = parts[5]
                cov = int(float(parts[9])) # ENCODE bedMethyl column 10
                pct = float(parts[10])     # ENCODE bedMethyl column 11
            except Exception:
                qc["parse_fail"] += 1
                continue

            qc["rows"] += 1
            if cov < 10:
                qc["coverage_lt10"] += 1
                continue

            seq = seqs.get(chrom)
            if seq is None:
                qc["chrom_missing"] += 1
                continue
            if start < 3 or start + 4 > len(seq):
                qc["edge"] += 1
                continue

            seven = seq[start-3:start+4]
            if strand == "-":
                seven = revcomp(seven)
            elif strand != "+":
                # Fallback for unknown strand: a reference G at the center
                # represents a cytosine on the opposite strand.
                if seven[3] == "G":
                    seven = revcomp(seven)

            if len(seven) != 7 or any(b not in BASES for b in seven):
                qc["non_acgt"] += 1
                continue
            if seven[3] != "C":
                qc["center_not_C"] += 1
                continue

            tri = seven[3:6]
            if context == "CHG":
                valid = (tri[1] != "G" and tri[2] == "G")
            elif context == "CHH":
                valid = (tri[1] != "G" and tri[2] != "G")
            else:
                raise ValueError("context must be CHG or CHH")
            if not valid:
                qc["context_mismatch"] += 1
                continue

            band = band_for(cov)
            qc["kept"] += 1

            # Overall >=10x plus four coverage strata.
            for bnd in ("ALL_10plus", band):
                add(overall, bnd, cov, pct)
                add(subtype, (bnd, tri), cov, pct)
                add(motif, (bnd, seven), cov, pct)
                for pos, idx in [(-3,0),(-2,1),(-1,2),(1,4),(2,5),(3,6)]:
                    add(position, (bnd, pos, seven[idx]), cov, pct)

            # 10%-wide methylation bins: 0,10,...,90,100.
            mbin = min(100, int(pct // 10) * 10)
            add(bin_subtype, (mbin, tri), cov, pct)
            add(bin_motif, (mbin, seven), cov, pct)
            for pos, idx in [(-3,0),(-2,1),(-1,2),(1,4),(2,5),(3,6)]:
                add(bin_position, (mbin, pos, seven[idx]), cov, pct)

    prefix = outdir / f"{label}_{context}"
    write_agg(str(prefix) + "_overall.csv",
              ["coverage_band"], overall, lambda k:(k,))
    write_agg(str(prefix) + "_subtypes.csv",
              ["coverage_band","subtype"], subtype, lambda k:k)
    write_agg(str(prefix) + "_position_summary.csv",
              ["coverage_band","position","base"], position, lambda k:k)
    write_agg(str(prefix) + "_7mer.csv",
              ["coverage_band","motif"], motif, lambda k:k)
    write_agg(str(prefix) + "_methylation_bin_subtypes.csv",
              ["methylation_bin_start_percent","subtype"], bin_subtype, lambda k:k)
    write_agg(str(prefix) + "_methylation_bin_positions.csv",
              ["methylation_bin_start_percent","position","base"], bin_position, lambda k:k)
    write_agg(str(prefix) + "_methylation_bin_7mer.csv",
              ["methylation_bin_start_percent","motif"], bin_motif, lambda k:k)

    qc_path = str(prefix) + "_QC.csv"
    with open(qc_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        for k in sorted(qc):
            w.writerow([k, qc[k]])
        w.writerow(["elapsed_minutes", f"{(time.time()-t0)/60:.3f}"])

    print(f"Done: {label} {context}. Kept {qc['kept']:,} sites.")
    print(f"Outputs: {prefix}_*.csv")

def main():
    print("Motor-neuron CpH sequence-context reducer v1")
    print("ENCODE ENCSR778ACS | GRCh38 | CHG + CHH | two isogenic replicates")

    ref = clean_path(input("\nGRCh38 primary-assembly FASTA (.fa.gz):\n> "))
    if not os.path.exists(ref):
        sys.exit(f"Reference not found: {ref}")

    default_dir = r"C:\IGOR\Papers\CG methylation contexts\Motor_neuron_WGBS"
    data_dir = clean_path(input(
        f"\nFolder containing the four downloaded BED files\n[{default_dir}]:\n> "
    )) or default_dir
    if not os.path.isdir(data_dir):
        sys.exit(f"Folder not found: {data_dir}")

    files = [
        ("ENCFF093TBN.bed.gz", "MotorNeuron_rep1", "CHG"),
        ("ENCFF109XWQ.bed.gz", "MotorNeuron_rep2", "CHG"),
        ("ENCFF505WSM.bed.gz", "MotorNeuron_rep1", "CHH"),
        ("ENCFF437XEW.bed.gz", "MotorNeuron_rep2", "CHH"),
    ]

    missing = [fn for fn,_,_ in files if not os.path.exists(os.path.join(data_dir, fn))]
    if missing:
        print("\nThese expected files were not found:")
        for fn in missing:
            print(" ", fn)
        print("\nCheck that the files are in the folder above and still end in .bed.gz.")
        sys.exit(1)

    outdir = Path(data_dir) / "Motor_neuron_CpH_context_results"
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = load_reference(ref)

    for fn, label, context in files:
        process_bed(os.path.join(data_dir, fn), label, context, seqs, outdir)

    print("\nALL FOUR FILES FINISHED.")
    print(f"Result folder:\n{outdir}")
    print("\nZIP that result folder and upload the ZIP to ChatGPT.")
    print("Do NOT upload the four large BED files.")

if __name__ == "__main__":
    main()
