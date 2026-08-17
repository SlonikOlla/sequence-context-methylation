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
            try: yield p[0],int(p[1]),p[5],float(p[9]),float(p[10])
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

class Stat:
    __slots__=("n","s","ss","cov")
    def __init__(self): self.n=0; self.s=0.; self.ss=0.; self.cov=0.
    def add(self,m,c): self.n+=1; self.s+=m; self.ss+=m*m; self.cov+=c
    def row(self):
        mean=self.s/self.n
        var=max(0.,(self.ss-self.s*self.s/self.n)/(self.n-1)) if self.n>1 else 0.
        return [self.n,mean,var**0.5,self.cov/self.n]

def covband(cov):
    if 10 <= cov < 15: return "10-14x"
    if 15 <= cov < 20: return "15-19x"
    if 20 <= cov < 30: return "20-29x"
    if cov >= 30: return "30plus"
    return None

def main():
    print("\nCHG coverage-stratified sequence-context reducer v1\n")
    genome=ask_file("GRCh38 primary-assembly FASTA (.fa.gz):\n> ")
    bed=ask_file("ENCODE CHG bedMethyl (.bed.gz):\n> ")
    sample=input("Short sample label (e.g. H1_CHG_rep1): ").strip() or "CHG_sample"
    sample=re.sub(r"[^A-Za-z0-9_-]+","_",sample)

    out=Path.home()/"Downloads"/f"WGBS_CHG_coverage_{sample}"
    out.mkdir(parents=True,exist_ok=True)
    t=out/"_write_test.tmp"; t.write_text("ok"); t.unlink()
    print("Output:",out)

    # Cache canonical chromosomes once; no preliminary pass through the huge BED.
    canonical={f"chr{i}" for i in range(1,23)}|{"chrX","chrY","chrM"}
    cache=Path(tempfile.mkdtemp(prefix="grch38_chgcov_"))
    print("\nBuilding chromosome cache...")
    found=set()
    for chrom,seq in fasta_records(genome):
        if chrom in canonical:
            with gzip.open(cache/f"{chrom}.seq.gz","wt",compresslevel=1) as f: f.write(seq)
            found.add(chrom); print(" ",chrom,len(seq))

    motif=defaultdict(Stat)     # band,motif
    pos=defaultdict(Stat)       # band,pos,base
    subtype=defaultdict(Stat)   # band,subtype
    overall=defaultdict(Stat)   # band
    qc=defaultdict(int)
    current=None; seq=None

    print("\nScanning CHG BED once...")
    try:
        for chrom,p,strand,cov,meth in bed_rows(bed):
            qc["total_rows"]+=1
            if chrom!=current:
                current=chrom
                cp=cache/f"{chrom}.seq.gz"
                if cp.exists():
                    with gzip.open(cp,"rt") as f: seq=f.read().strip().upper()
                else: seq=None
                print(" ",chrom)
            if seq is None:
                qc["missing_chrom"]+=1; continue
            band=covband(cov)
            if band is None:
                qc["below10"]+=1; continue
            got=context_for(seq,p,strand)
            if got is None:
                qc["bad_context"]+=1; continue
            ctx,sub,vars_=got
            qc["kept"]+=1
            for cb in (band,"ALL_10plus"):
                overall[cb].add(meth,cov)
                motif[(cb,ctx)].add(meth,cov)
                subtype[(cb,sub)].add(meth,cov)
                for pp,b in vars_: pos[(cb,pp,b)].add(meth,cov)
    finally:
        shutil.rmtree(cache,ignore_errors=True)

    with open(out/"QC.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["metric","value"])
        for k in ["total_rows","kept","below10","bad_context","missing_chrom"]: w.writerow([k,qc[k]])

    def write_stats(fn,header,items):
        with open(out/fn,"w",newline="") as f:
            w=csv.writer(f); w.writerow(header)
            for keys,st in items:
                if not isinstance(keys,tuple): keys=(keys,)
                w.writerow(list(keys)+st.row())

    write_stats("coverage_overall.csv",
        ["coverage_band","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(overall.items()))
    write_stats("coverage_CHG_subtypes.csv",
        ["coverage_band","subtype","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(subtype.items()))
    write_stats("coverage_position_summary.csv",
        ["coverage_band","position","base","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(pos.items()))
    write_stats("coverage_motif_7mer.csv",
        ["coverage_band","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(motif.items()))

    selected={"TTACAGC","CTACAGC","TTACAGG","TTACAGT","TTACAGA"}
    with open(out/"coverage_selected_motifs.csv","w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["coverage_band","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
        for (cb,m),st in sorted(motif.items()):
            if m in selected: w.writerow([cb,m]+st.row())

    (out/"README.txt").write_text(
        "CHG coverage-stratified sequence-context analysis\n"
        "Coverage strata: 10-14x, 15-19x, 20-29x, >=30x, plus ALL_10plus.\n"
        "ZIP this folder and upload it to ChatGPT.\n",
        encoding="utf-8"
    )
    print("\nDONE")
    print("Kept CHG sites >=10x:",f"{qc['kept']:,}")
    print("Results:",out)

if __name__=="__main__":
    main()
