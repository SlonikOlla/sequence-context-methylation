#!/usr/bin/env python3
"""
Robust reducer for ENCODE WGBS bedMethyl + GRCh38 sequence context.

Fix in v2:
ENCODE BED files can be ordered chr1, chr10, chr11, ... while GENCODE FASTA is
chr1, chr2, chr3, ... . This version does not assume matching chromosome order.

It first makes a temporary compressed chromosome cache from the FASTA, then
streams each BED file chromosome-by-chromosome. Only one chromosome sequence is
held in memory at a time.

Outputs:
  *_QC.csv
  *_motif_summary.csv
  *_position_summary.csv
  *_bin_position_counts.csv
  *_chromosome_summary.csv
"""
import csv, gzip, os, shutil, sys, tempfile
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = set("ACGT")

def revcomp(s):
    return s.translate(COMP)[::-1].upper()

def op_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).lower().endswith(".gz") else open(path, mode)

def fasta_records(path):
    with op_text(path, "rt") as f:
        name=None; seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name=line[1:].split()[0]
                seq=[]
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq).upper()

def bed_rows(path):
    with op_text(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            p=line.rstrip("\n").split("\t")
            if len(p)<11:
                p=line.rstrip("\n").split()
            if len(p)<11:
                continue
            try:
                chrom=p[0]
                start=int(p[1])
                strand=p[5]
                cov=float(p[9])
                meth=float(p[10])
            except Exception:
                continue
            yield chrom,start,strand,cov,meth

def get_bed_chromosomes(paths):
    chroms=set()
    for path in paths:
        print(f"Reading chromosome names from {Path(path).name} ...")
        for chrom,_,_,_,_ in bed_rows(path):
            chroms.add(chrom)
    return chroms

def build_chrom_cache(genome, needed, cache_dir):
    print("\nBuilding temporary chromosome cache from GRCh38...")
    found=set()
    for chrom,seq in fasta_records(genome):
        if chrom not in needed:
            continue
        out=cache_dir/f"{chrom.replace('/','_')}.seq.gz"
        with gzip.open(out,"wt",compresslevel=1) as f:
            f.write(seq)
        found.add(chrom)
        print(f"  cached {chrom} ({len(seq):,} bp)")
    missing=needed-found
    if missing:
        print("\nWARNING: BED chromosomes absent from FASTA:")
        print(", ".join(sorted(missing)[:50]))
    return found

def load_chrom(cache_dir, chrom):
    p=cache_dir/f"{chrom.replace('/','_')}.seq.gz"
    if not p.exists():
        return None
    with gzip.open(p,"rt") as f:
        return f.read().strip().upper()

def context_for(seq,p,strand,label):
    L=len(seq)
    if label=="CpG":
        if strand=="+":
            if p<3 or p+5>L: return None
            s=seq[p-3:p+5]
        elif strand=="-":
            if p<4 or p+4>L: return None
            s=revcomp(seq[p-4:p+4])
        else:
            # Some bedMethyl outputs may use ".". Try plus orientation first,
            # then minus if the reference base indicates it.
            if p>=3 and p+5<=L:
                s=seq[p-3:p+5]
                if len(s)==8 and s[3:5]=="CG":
                    pass
                else:
                    if p<4 or p+4>L: return None
                    s=revcomp(seq[p-4:p+4])
            else:
                return None
        if len(s)!=8 or s[3:5]!="CG" or any(c not in BASES for c in s):
            return None
        vars_=[("-3",s[0]),("-2",s[1]),("-1",s[2]),("+1",s[5]),("+2",s[6]),("+3",s[7])]
        return s,vars_

    if strand=="+":
        if p<3 or p+4>L: return None
        s=seq[p-3:p+4]
    elif strand=="-":
        if p<3 or p+4>L: return None
        s=revcomp(seq[p-3:p+4])
    else:
        if p<3 or p+4>L: return None
        s=seq[p-3:p+4]
        if s[3]!="C":
            s=revcomp(seq[p-3:p+4])

    if len(s)!=7 or s[3]!="C" or any(c not in BASES for c in s):
        return None
    if label=="CHG":
        if s[4]=="G" or s[5]!="G": return None
    elif label=="CHH":
        if s[4]=="G" or s[5]=="G": return None
    vars_=[("-3",s[0]),("-2",s[1]),("-1",s[2]),("+1",s[4]),("+2",s[5]),("+3",s[6])]
    return s,vars_

class Agg:
    def __init__(self,label,mincov):
        self.label=label; self.mincov=mincov
        self.total=0; self.kept=0; self.low=0; self.bad=0; self.missing_chrom=0
        self.motif=defaultdict(lambda:[0,0.,0.,0.])
        self.pos=defaultdict(lambda:[0,0.,0.,0.])
        self.binpos=defaultdict(int)
        self.chrom=defaultdict(lambda:[0,0.,0.,0.])

    def add(self,chrom,context,meth,cov,vars_):
        self.kept+=1
        x=self.motif[context]; x[0]+=1; x[1]+=meth; x[2]+=meth*meth; x[3]+=cov
        x=self.chrom[chrom]; x[0]+=1; x[1]+=meth; x[2]+=meth*meth; x[3]+=cov
        b=min(9,int(meth//10))
        binlab=f"{b*10}-{100 if b==9 else (b+1)*10}%"
        for pos,base in vars_:
            x=self.pos[(pos,base)]; x[0]+=1; x[1]+=meth; x[2]+=meth*meth; x[3]+=cov
            self.binpos[(binlab,pos,base)]+=1

def process_bed(path,label,mincov,cache_dir,found):
    agg=Agg(label,mincov)
    current=None; seq=None
    print(f"\nProcessing {Path(path).name} as {label}...")
    for chrom,p,strand,cov,meth in bed_rows(path):
        agg.total+=1
        if chrom!=current:
            current=chrom
            if chrom in found:
                seq=load_chrom(cache_dir,chrom)
                print(f"  {chrom}")
            else:
                seq=None
        if seq is None:
            agg.missing_chrom+=1
            continue
        if cov<mincov:
            agg.low+=1
            continue
        got=context_for(seq,p,strand,label)
        if got is None:
            agg.bad+=1
            continue
        context,vars_=got
        agg.add(chrom,context,meth,cov,vars_)
    return agg

def stats(v):
    n,s,ss,sc=v
    mean=s/n
    var=max(0.,(ss-s*s/n)/(n-1)) if n>1 else 0.
    return n,mean,var**0.5,sc/n

def cleanstem(path):
    stem=Path(path).name
    for suf in [".gz",".bed",".bedMethyl",".bedmethyl"]:
        if stem.endswith(suf):
            stem=stem[:-len(suf)]
    return stem

def write_outputs(agg,path,outdir):
    pre=outdir/cleanstem(path)
    with open(str(pre)+"_QC.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["metric","value"])
        for a,b in [
            ("context_class",agg.label),("coverage_threshold",agg.mincov),
            ("total_bed_rows_seen",agg.total),("rows_kept",agg.kept),
            ("rows_below_coverage",agg.low),
            ("rows_bad_or_mismatched_context",agg.bad),
            ("rows_on_chromosomes_absent_from_fasta",agg.missing_chrom)]:
            w.writerow([a,b])

    with open(str(pre)+"_motif_summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["context_class","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
        for motif,v in sorted(agg.motif.items()):
            n,m,sd,c=stats(v); w.writerow([agg.label,motif,n,m,sd,c])

    orderpos={"-3":0,"-2":1,"-1":2,"+1":3,"+2":4,"+3":5}
    with open(str(pre)+"_position_summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["context_class","position","base","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
        for (pos,base),v in sorted(agg.pos.items(),key=lambda x:(orderpos[x[0][0]],x[0][1])):
            n,m,sd,c=stats(v); w.writerow([agg.label,pos,base,n,m,sd,c])

    orderbin={f"{i*10}-{(i+1)*10}%":i for i in range(9)}; orderbin["90-100%"]=9
    with open(str(pre)+"_bin_position_counts.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["context_class","methylation_bin","position","base","site_count"])
        for (b,pos,base),n in sorted(agg.binpos.items(),key=lambda x:(orderbin[x[0][0]],orderpos[x[0][1]],x[0][2])):
            w.writerow([agg.label,b,pos,base,n])

    with open(str(pre)+"_chromosome_summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["context_class","chromosome","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
        for chrom,v in agg.chrom.items():
            n,m,sd,c=stats(v); w.writerow([agg.label,chrom,n,m,sd,c])

def ask_file(msg):
    while True:
        p=input(msg).strip().strip('"')
        if os.path.isfile(p): return p
        print("File not found. Try again.")

def main():
    print("\nWGBS sequence-context reducer v3\n")
    genome=ask_file("Paste full path to GRCh38 primary-assembly FASTA (.fa.gz):\n> ")
    try: mincov=int(input("Minimum coverage [10]: ").strip() or "10")
    except: mincov=10

    entries=[]
    print("\nAdd bedMethyl files. Press Enter at file prompt when finished.")
    while True:
        p=input("\nbedMethyl file (.bed.gz), or Enter to finish:\n> ").strip().strip('"')
        if not p: break
        if not os.path.isfile(p):
            print("File not found."); continue
        lab=input("Context class [CpG/CHG/CHH]: ").strip()
        lab={"cpg":"CpG","chg":"CHG","chh":"CHH"}.get(lab.lower())
        if not lab:
            print("Invalid context class."); continue
        entries.append((p,lab))

    if not entries:
        print("No BED files supplied."); return

    outdir=Path.home()/"Downloads"/"WGBS_context_results_v3"
    outdir.mkdir(parents=True, exist_ok=True)
    # Fail early if Windows cannot write here, before the expensive scan.
    testfile=outdir/"_write_test.tmp"
    try:
        testfile.write_text("ok", encoding="utf-8")
        testfile.unlink()
    except Exception as e:
        print("\nCannot write to output folder:", outdir)
        print("Error:", e)
        return
    print("\nOutput folder:", outdir)
    cache=Path(tempfile.mkdtemp(prefix="grch38_chrom_cache_"))

    try:
        needed=get_bed_chromosomes([p for p,_ in entries])
        found=build_chrom_cache(genome,needed,cache)
        for path,label in entries:
            agg=process_bed(path,label,mincov,cache,found)
            write_outputs(agg,path,outdir)
            print(f"\n{Path(path).name}:")
            print(f"  total BED rows: {agg.total:,}")
            print(f"  kept at >= {mincov}x: {agg.kept:,}")
            print(f"  below coverage: {agg.low:,}")
            print(f"  context mismatches: {agg.bad:,}")
            print(f"  absent chromosome rows: {agg.missing_chrom:,}")
    finally:
        shutil.rmtree(cache,ignore_errors=True)

    with open(outdir/"README.txt","w") as f:
        f.write(f"Coverage threshold: {mincov}\n")
        f.write("Robust chromosome-order-independent reducer v3.\n")
        f.write("ZIP this folder and upload it to ChatGPT.\n")

    print("\nDone.")
    print("Results folder:")
    print(outdir)
    print("\nZIP that folder and upload it here.")

if __name__=="__main__":
    main()
