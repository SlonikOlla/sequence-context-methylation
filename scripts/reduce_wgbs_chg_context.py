#!/usr/bin/env python3
import csv, gzip, os, re, shutil, tempfile
from pathlib import Path
from collections import defaultdict

BASES=set("ACGT")
COMP=str.maketrans("ACGTNacgtn","TGCANtgcan")

def revcomp(s): return s.translate(COMP)[::-1].upper()
def op_text(path,mode="rt"):
    return gzip.open(path,mode) if str(path).lower().endswith(".gz") else open(path,mode)

def fasta_records(path):
    with op_text(path,"rt") as f:
        name=None; seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None: yield name,"".join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line.strip())
        if name is not None: yield name,"".join(seq).upper()

def bed_rows(path):
    with op_text(path,"rt") as f:
        for line in f:
            if not line.strip() or line.startswith(("#","track","browser")): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<11: p=line.rstrip("\n").split()
            if len(p)<11: continue
            try:
                yield p[0],int(p[1]),p[5],float(p[9]),float(p[10])
            except: continue

def ask_file(msg):
    while True:
        p=input(msg).strip().strip('"')
        if os.path.isfile(p): return p
        print("File not found. Try again.")

def context_for(seq,p,strand):
    L=len(seq)
    if p<3 or p+4>L: return None
    s=seq[p-3:p+4]
    if strand=="-": s=revcomp(s)
    elif strand not in ("+","."): return None
    elif strand=="." and s[3]!="C": s=revcomp(s)
    if len(s)!=7 or s[3]!="C" or any(c not in BASES for c in s): return None
    if s[4]=="G" or s[5]!="G": return None
    subtype=s[3:6]
    vars_=[("-3",s[0]),("-2",s[1]),("-1",s[2]),("+1",s[4]),("+2",s[5]),("+3",s[6])]
    return s,subtype,vars_

def stats(v):
    n,s,ss,sc=v
    mean=s/n
    var=max(0.,(ss-s*s/n)/(n-1)) if n>1 else 0.
    return n,mean,var**0.5,sc/n

def main():
    print("\nCHG sequence-context reducer v1\n")
    genome=ask_file("GRCh38 primary-assembly FASTA (.fa.gz):\n> ")
    bed=ask_file("ENCODE CHG bedMethyl (.bed.gz):\n> ")
    try: mincov=int(input("Minimum coverage [10]: ").strip() or "10")
    except: mincov=10
    sample=input("Short sample label [H1_CHG]: ").strip() or "H1_CHG"
    sample=re.sub(r"[^A-Za-z0-9_-]+","_",sample)

    outdir=Path.home()/"Downloads"/f"WGBS_CHG_results_{sample}"
    outdir.mkdir(parents=True,exist_ok=True)
    test=outdir/"_write_test.tmp"; test.write_text("ok"); test.unlink()

    needed={c for c,_,_,_,_ in bed_rows(bed)}
    cache=Path(tempfile.mkdtemp(prefix="grch38_chg_cache_"))
    found=set()
    print("\nBuilding temporary chromosome cache...")
    for chrom,seq in fasta_records(genome):
        if chrom in needed:
            with gzip.open(cache/f"{chrom}.seq.gz","wt",compresslevel=1) as f: f.write(seq)
            found.add(chrom); print(f"  cached {chrom} ({len(seq):,} bp)")

    motif=defaultdict(lambda:[0,0.,0.,0.]); pos=defaultdict(lambda:[0,0.,0.,0.])
    bins=defaultdict(int); subtype=defaultdict(lambda:[0,0.,0.,0.]); chromsum=defaultdict(lambda:[0,0.,0.,0.])
    total=kept=low=bad=missing=0; current=None; seq=None

    try:
        print("\nProcessing CHG bedMethyl...")
        for chrom,p,strand,cov,meth in bed_rows(bed):
            total+=1
            if chrom!=current:
                current=chrom
                cp=cache/f"{chrom}.seq.gz"
                if cp.exists():
                    with gzip.open(cp,"rt") as f: seq=f.read().strip().upper()
                else: seq=None
                print(" ",chrom)
            if seq is None: missing+=1; continue
            if cov<mincov: low+=1; continue
            got=context_for(seq,p,strand)
            if got is None: bad+=1; continue
            ctx,sub,vars_=got; kept+=1
            for d,k in [(motif,ctx),(subtype,sub),(chromsum,chrom)]:
                x=d[k]; x[0]+=1; x[1]+=meth; x[2]+=meth*meth; x[3]+=cov
            b=min(9,int(meth//10)); binlab=f"{b*10}-{100 if b==9 else (b+1)*10}%"
            for pp,base in vars_:
                x=pos[(pp,base)]; x[0]+=1; x[1]+=meth; x[2]+=meth*meth; x[3]+=cov
                bins[(binlab,pp,base)]+=1

        with open(outdir/"QC.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["metric","value"])
            for a,b in [("coverage_threshold",mincov),("total_rows",total),("kept",kept),("below_coverage",low),("bad_context",bad),("missing_chrom",missing)]: w.writerow([a,b])

        with open(outdir/"CHG_subtypes.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["subtype","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
            for k,v in sorted(subtype.items()): w.writerow([k,*stats(v)])

        order={"-3":0,"-2":1,"-1":2,"+1":3,"+2":4,"+3":5}
        with open(outdir/"position_summary.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["position","base","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
            for (pp,b),v in sorted(pos.items(),key=lambda x:(order[x[0][0]],x[0][1])): w.writerow([pp,b,*stats(v)])

        with open(outdir/"motif_7mer.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
            for k,v in sorted(motif.items()): w.writerow([k,*stats(v)])

        orderbin={f"{i*10}-{(i+1)*10}%":i for i in range(9)}; orderbin["90-100%"]=9
        with open(outdir/"bin_position_counts.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["methylation_bin","position","base","site_count"])
            for (bb,pp,b),n in sorted(bins.items(),key=lambda x:(orderbin[x[0][0]],order[x[0][1]],x[0][2])): w.writerow([bb,pp,b,n])

        with open(outdir/"chromosome_summary.csv","w",newline="") as f:
            w=csv.writer(f); w.writerow(["chromosome","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
            for k,v in chromsum.items(): w.writerow([k,*stats(v)])

        (outdir/"README.txt").write_text(f"CHG sequence-context reducer v1\nSample: {sample}\nCoverage threshold: {mincov}\n")
    finally:
        shutil.rmtree(cache,ignore_errors=True)

    print("\nDone.\nResults folder:\n",outdir,"\nZIP this folder and upload it here.")

if __name__=="__main__":
    main()
