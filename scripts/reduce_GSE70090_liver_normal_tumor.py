#!/usr/bin/env python3
"""
GSE70090 matched liver normal-tumor CpG sequence-context reducer.

Inputs
------
1) hg19.fa.gz
2) GSE70090_RAW.tar

The TAR contains 28 GEO processed CpG files. This script automatically selects
the 16 liver files:
  N1-N4 and T1-T4 x BS and oxBS.

For each file it:
  * reads methylated-read count and total-read count at each CpG
  * detects whether the CpG coordinate is 0- or 1-based by checking hg19
  * requires coverage >=10
  * extracts the strand-normalized 8-mer NNNCGNNN and 6-mer NNCGNN
  * summarizes methylation by:
      - complete 6-mer
      - +/-3 base position
      - coverage band
      - chromosome
  * writes a compact site sample for later hierarchical modeling

It also writes matched normal-vs-tumor summaries for BS and oxBS and
aggregate BS-minus-oxBS estimates (a proxy for 5hmC contribution).

No external Python packages are required.
"""

import csv, gzip, io, os, re, sys, tarfile, time, random
from collections import defaultdict
from pathlib import Path

BASES=set("ACGT")
COMP=str.maketrans("ACGTNacgtn","TGCANtgcan")

def revcomp(s):
    return s.translate(COMP)[::-1].upper()

def clean(s):
    return s.strip().strip('"').strip("'")

def open_text(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path,"rt",encoding="ascii",errors="ignore")
    return open(path,"rt",encoding="ascii",errors="ignore")

def fasta(path):
    with open_text(path) as f:
        name=None; seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name,"".join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            else:
                seq.append(line.strip())
        if name is not None:
            yield name,"".join(seq).upper()

def load_reference(path):
    print("\nLoading hg19 reference...")
    d={}
    for chrom,seq in fasta(path):
        d[chrom]=seq
        if chrom in {f"chr{i}" for i in range(1,23)}|{"chrX","chrY","chrM"}:
            print(f"  {chrom}: {len(seq):,} bp")
    print(f"Loaded {len(d)} reference sequences.")
    return d

def band(cov):
    if cov < 10: return None
    if cov <=14: return "10-14x"
    if cov <=19: return "15-19x"
    if cov <=29: return "20-29x"
    return "30plus"

class Stat:
    __slots__=("n","sum_pct","sum_cov","sum_meth")
    def __init__(self):
        self.n=0; self.sum_pct=0.0; self.sum_cov=0; self.sum_meth=0
    def add(self,m,c):
        self.n+=1
        self.sum_cov+=c
        self.sum_meth+=m
        self.sum_pct+=100.0*m/c
    def vals(self):
        return (
            self.n,
            self.sum_cov,
            self.sum_meth,
            self.sum_pct/self.n if self.n else float("nan"),
            100.0*self.sum_meth/self.sum_cov if self.sum_cov else float("nan"),
        )

def classify_name(name):
    # Example: GSM1716957_liver_N1_BS.cg.txt.gz
    b=os.path.basename(name)
    m=re.search(r"liver_([NT])([1-4])_(oxBS|BS)\.cg\.txt\.gz$",b,re.I)
    if not m: return None
    status={"N":"Normal","T":"Tumor"}[m.group(1).upper()]
    patient=int(m.group(2))
    assay="oxBS" if m.group(3).lower()=="oxbs" else "BS"
    return status,patient,assay,b

def parse_numeric_row(parts):
    # GEO states: methylated read count and total read count at each CpG.
    # Most files are expected to have chromosome, coordinate and then count fields.
    if len(parts)<4: return None
    chrom=parts[0]
    if not chrom.startswith("chr"): return None
    try:
        pos=int(float(parts[1]))
    except:
        return None

    nums=[]
    for idx in range(2,len(parts)):
        try:
            v=int(float(parts[idx]))
            nums.append((idx,v))
        except:
            continue

    # choose a pair m,total among numeric columns with 0<=m<=total.
    # Prefer adjacent numeric columns nearest the end.
    candidates=[]
    for i in range(len(nums)):
        for j in range(i+1,len(nums)):
            mi,m=nums[i]; ti,t=nums[j]
            if 0 <= m <= t and t>=0:
                score=(ti-mi==1, ti, -mi)
                candidates.append((score,m,t))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    _,m,t=candidates[0]
    return chrom,pos,m,t

def detect_coordinate_mode(tar,member,seqs):
    # Determine whether GEO coordinate points to C as 0-based or 1-based.
    fh=tar.extractfile(member)
    gz=gzip.GzipFile(fileobj=fh,mode="rb")
    txt=io.TextIOWrapper(gz,encoding="ascii",errors="ignore")
    scores={"zero":0,"one":0}
    n=0
    for line in txt:
        p=line.strip().split()
        row=parse_numeric_row(p)
        if not row: continue
        chrom,pos,m,t=row
        seq=seqs.get(chrom)
        if seq is None: continue
        if 0 <= pos < len(seq)-1 and seq[pos:pos+2]=="CG":
            scores["zero"]+=1
        if 1 <= pos <= len(seq)-1 and seq[pos-1:pos+1]=="CG":
            scores["one"]+=1
        n+=1
        if n>=3000: break
    mode=max(scores,key=scores.get)
    print(f"  coordinate detection: {scores} -> {mode}-based")
    if scores[mode] < max(20, int(0.5*n)):
        print("  WARNING: coordinate matching is weaker than expected.")
    return mode

def context(seq,pos,mode):
    c0=pos if mode=="zero" else pos-1
    if c0<3 or c0+5>len(seq): return None
    s=seq[c0-3:c0+5].upper()
    if len(s)!=8 or s[3:5]!="CG" or any(x not in BASES for x in s):
        return None
    six=s[1:7]
    positions=[(-3,s[0]),(-2,s[1]),(-1,s[2]),(1,s[5]),(2,s[6]),(3,s[7])]
    return s,six,positions

def write_stats(path,header,agg):
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(header+["site_count","total_coverage","methylated_reads",
                           "mean_methylation_percent","coverage_weighted_methylation_percent"])
        for k in sorted(agg,key=lambda z:str(z)):
            keys=k if isinstance(k,tuple) else (k,)
            w.writerow(list(keys)+list(agg[k].vals()))

def process_member(tar,member,meta,seqs,outdir,reservoir_n=250000):
    status,patient,assay,bname=meta
    label=f"P{patient}_{status}_{assay}"
    print(f"\n=== {label} ===")
    mode=detect_coordinate_mode(tar,member,seqs)

    overall=defaultdict(Stat)
    sixmer=defaultdict(Stat)
    posagg=defaultdict(Stat)
    chromagg=defaultdict(Stat)
    qc=defaultdict(int)

    # Uniform reservoir sample for later site-level modeling.
    rng=random.Random(1000+patient+(0 if status=="Normal" else 100)+(0 if assay=="BS" else 200))
    reservoir=[]
    seen_kept=0

    fh=tar.extractfile(member)
    gz=gzip.GzipFile(fileobj=fh,mode="rb")
    txt=io.TextIOWrapper(gz,encoding="ascii",errors="ignore")

    t0=time.time()
    for i,line in enumerate(txt,1):
        if i%5_000_000==0:
            print(f"  {i:,} rows; {qc['kept']:,} kept; {(time.time()-t0)/60:.1f} min")
        parts=line.strip().split()
        row=parse_numeric_row(parts)
        if not row:
            qc["parse_or_header"]+=1; continue
        chrom,p,meth,cov=row
        qc["rows"]+=1
        if cov<10:
            qc["coverage_lt10"]+=1; continue
        seq=seqs.get(chrom)
        if seq is None:
            qc["chrom_missing"]+=1; continue
        got=context(seq,p,mode)
        if got is None:
            qc["context_mismatch"]+=1; continue
        eight,six,positions=got
        qc["kept"]+=1
        bnd=band(cov)
        pct=100.0*meth/cov

        for bb in ("ALL_10plus",bnd):
            overall[bb].add(meth,cov)
            sixmer[(bb,six)].add(meth,cov)
            for pp,base in positions:
                posagg[(bb,pp,base)].add(meth,cov)
        chromagg[chrom].add(meth,cov)

        seen_kept+=1
        rec=(chrom,p,meth,cov,pct,bnd,six,eight,
             eight[0],eight[1],eight[2],eight[5],eight[6],eight[7])
        if len(reservoir)<reservoir_n:
            reservoir.append(rec)
        else:
            j=rng.randrange(seen_kept)
            if j<reservoir_n:
                reservoir[j]=rec

    prefix=outdir/label
    write_stats(str(prefix)+"_overall.csv",["coverage_band"],overall)
    write_stats(str(prefix)+"_6mer.csv",["coverage_band","motif"],sixmer)
    write_stats(str(prefix)+"_position.csv",["coverage_band","position","base"],posagg)
    write_stats(str(prefix)+"_chromosome.csv",["chromosome"],chromagg)

    with open(str(prefix)+"_site_sample.csv","w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["patient","status","assay","chrom","coordinate_raw","methylated_reads","coverage",
                    "methylation_percent","coverage_band","sixmer","eightmer",
                    "minus3","minus2","minus1","plus1","plus2","plus3"])
        for r in reservoir:
            w.writerow([patient,status,assay]+list(r))

    with open(str(prefix)+"_QC.csv","w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["metric","value"])
        w.writerow(["source_member",member.name])
        w.writerow(["coordinate_mode",mode])
        for k in sorted(qc): w.writerow([k,qc[k]])
        w.writerow(["reservoir_sample_size",len(reservoir)])
        w.writerow(["elapsed_minutes",f"{(time.time()-t0)/60:.3f}"])

    print(f"Finished {label}: {qc['kept']:,} CpGs >=10x")
    return label

def read_6mer(path):
    d={}
    with open(path,newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["coverage_band"]!="ALL_10plus": continue
            d[r["motif"]]={
                "n":int(r["site_count"]),
                "mean":float(r["mean_methylation_percent"]),
                "weighted":float(r["coverage_weighted_methylation_percent"])
            }
    return d

def make_pair_summaries(outdir):
    # Normal vs Tumor within patient and assay, plus BS-oxBS at aggregate motif level.
    rows=[]
    for patient in range(1,5):
        for assay in ["BS","oxBS"]:
            n=read_6mer(outdir/f"P{patient}_Normal_{assay}_6mer.csv")
            t=read_6mer(outdir/f"P{patient}_Tumor_{assay}_6mer.csv")
            for motif in sorted(set(n)&set(t)):
                rows.append([patient,assay,motif,n[motif]["n"],t[motif]["n"],
                             n[motif]["mean"],t[motif]["mean"],t[motif]["mean"]-n[motif]["mean"]])
    with open(outdir/"paired_normal_tumor_6mer.csv","w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["patient","assay","motif","normal_sites","tumor_sites",
                    "normal_mean_methylation","tumor_mean_methylation","tumor_minus_normal_pp"])
        w.writerows(rows)

    rows=[]
    for patient in range(1,5):
        for status in ["Normal","Tumor"]:
            bs=read_6mer(outdir/f"P{patient}_{status}_BS_6mer.csv")
            ox=read_6mer(outdir/f"P{patient}_{status}_oxBS_6mer.csv")
            for motif in sorted(set(bs)&set(ox)):
                # BS approx 5mC+5hmC; oxBS approx 5mC.
                hm=bs[motif]["mean"]-ox[motif]["mean"]
                rows.append([patient,status,motif,bs[motif]["n"],ox[motif]["n"],
                             bs[motif]["mean"],ox[motif]["mean"],hm])
    with open(outdir/"BS_minus_oxBS_6mer_proxy5hmC.csv","w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["patient","status","motif","BS_sites","oxBS_sites",
                    "BS_mean_percent","oxBS_mean_percent","BS_minus_oxBS_pp"])
        w.writerows(rows)

def main():
    print("\nGSE70090 matched liver normal-tumor sequence-context reducer v1")
    print("4 matched liver pairs | BS + oxBS | hg19\n")

    ref=clean(input("hg19.fa.gz:\n> "))
    tarpath=clean(input("\nGSE70090_RAW.tar:\n> "))
    if not os.path.isfile(ref): sys.exit("hg19 reference not found.")
    if not os.path.isfile(tarpath): sys.exit("GSE70090_RAW.tar not found.")

    outdir=Path(tarpath).resolve().parent/"GSE70090_liver_context_results"
    outdir.mkdir(parents=True,exist_ok=True)
    t=outdir/"_write_test.tmp"; t.write_text("ok"); t.unlink()

    seqs=load_reference(ref)

    with tarfile.open(tarpath,"r") as tar:
        members=[]
        for m in tar.getmembers():
            if not m.isfile(): continue
            meta=classify_name(m.name)
            if meta:
                members.append((m,meta))
        members.sort(key=lambda z:(z[1][1],z[1][0],z[1][2]))

        print("\nDetected liver files:")
        for m,meta in members:
            print(" ",m.name)
        if len(members)!=16:
            print(f"\nWARNING: expected 16 liver files, detected {len(members)}.")
            print("If filenames differ, stop and send the list above to ChatGPT.")
            if len(members)==0: return

        for m,meta in members:
            process_member(tar,m,meta,seqs,outdir)

    make_pair_summaries(outdir)

    (outdir/"README.txt").write_text(
        "GSE70090 liver normal-tumor sequence-context results.\n"
        "BS measures bisulfite-resistant 5mC+5hmC; oxBS is used as an approximation of 5mC.\n"
        "BS-minus-oxBS motif summaries are exploratory proxies for 5hmC and are not site-matched estimates.\n"
        "Primary coverage threshold: >=10x.\n"
        "Upload this entire folder as a ZIP to ChatGPT.\n",
        encoding="utf-8"
    )

    print("\nALL 16 LIVER FILES FINISHED.")
    print("Results folder:")
    print(outdir)
    print("\nZIP that folder and upload it here.")

if __name__=="__main__":
    main()
