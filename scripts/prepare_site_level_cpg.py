#!/usr/bin/env python3
"""Prepare a balanced site-level CpG dataset for multivariable modeling.

Uses ENCODE bedMethyl, GRCh38, GENCODE GTF and UCSC CpG islands. Reservoir-
samples each SAMPLE x COVERAGE-BAND x PRIMARY-COMPARTMENT stratum, preserving
stratum population sizes so inverse-probability weights can be used downstream.
No third-party packages are required.
"""
import csv,gzip,os,re,tempfile,shutil,bisect,random
from pathlib import Path
from collections import defaultdict

BASES=set("ACGT"); COMP=str.maketrans("ACGTNacgtn","TGCANtgcan")
def revcomp(s):return s.translate(COMP)[::-1].upper()
def op_text(p,mode="rt"):
    return gzip.open(p,mode) if str(p).lower().endswith(".gz") else open(p,mode,encoding="utf-8")
def fasta_records(path):
    with op_text(path) as f:
        name=None;seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None:yield name,"".join(seq).upper()
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
        if name is not None:yield name,"".join(seq).upper()
def merge_intervals(iv):
    if not iv:return []
    iv=sorted(iv);out=[list(iv[0])]
    for s,e in iv[1:]:
        if s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
        else:out.append([s,e])
    return [(s,e) for s,e in out]
def point_in(m,starts,p):
    if not m:return False
    i=bisect.bisect_right(starts,p)-1
    return i>=0 and p<m[i][1]
def parse_attrs(s):return {m.group(1):m.group(2) for m in re.finditer(r'(\S+)\s+"([^"]*)"',s)}
def load_gtf(path):
    raw=defaultdict(lambda:defaultdict(list));print("Reading GTF...")
    with op_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):continue
            p=line.rstrip("\n").split("\t")
            if len(p)<9:continue
            chrom,_,feat,start,end,_,strand,_,attrs=p
            if not chrom.startswith("chr"):continue
            s=int(start)-1;e=int(end);a=parse_attrs(attrs);gt=a.get("gene_type") or a.get("gene_biotype") or ""
            if feat=="gene":
                raw[chrom]["gene"].append((s,e))
                tss=s if strand=="+" else e-1
                ps=max(0,tss-2000) if strand=="+" else max(0,tss-500)
                pe=tss+501 if strand=="+" else tss+2001
                raw[chrom]["promoter"].append((ps,pe))
            elif feat=="exon":raw[chrom]["exon"].append((s,e))
            elif feat=="CDS":raw[chrom]["CDS"].append((s,e))
    idx={}
    for c,types in raw.items():
        idx[c]={}
        for typ,iv in types.items():
            m=merge_intervals(iv);idx[c][typ]=(m,[x[0] for x in m])
    return idx
def load_cgi(path):
    raw=defaultdict(list);print("Reading CpG islands...")
    with op_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):continue
            p=line.rstrip("\n").split("\t")
            try:
                if p[0].isdigit() and p[1].startswith("chr"):c=p[1];s=int(p[2]);e=int(p[3])
                elif p[0].startswith("chr"):c=p[0];s=int(p[1]);e=int(p[2])
                else:continue
            except:continue
            raw[c].append((s,e))
    idx={}
    for c,iv in raw.items():
        m=merge_intervals(iv);idx[c]=(m,[x[0] for x in m])
    return idx
def flag(idx,c,t,p):
    x=idx.get(c,{}).get(t);return point_in(x[0],x[1],p) if x else False
def cgiflag(idx,c,p):
    x=idx.get(c);return point_in(x[0],x[1],p) if x else False
def primary(gidx,cidx,c,p):
    if cgiflag(cidx,c,p):return "CpG_island"
    if flag(gidx,c,"promoter",p):return "promoter_non_CGI"
    if flag(gidx,c,"CDS",p):return "CDS_nonpromoter_nonCGI"
    if flag(gidx,c,"exon",p):return "exon_noncoding_nonpromoter_nonCGI"
    if flag(gidx,c,"gene",p):return "intron_or_other_genic_nonpromoter_nonCGI"
    return "intergenic_nonCGI"
def bed_rows(path):
    with op_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#","track","browser")):continue
            p=line.rstrip("\n").split("\t")
            if len(p)<11:p=line.rstrip("\n").split()
            if len(p)<11:continue
            try:yield p[0],int(p[1]),p[5],int(round(float(p[9]))),float(p[10])
            except:continue
def cpg_context(seq,p,strand):
    L=len(seq)
    if strand=="+":
        if p<3 or p+5>L:return None
        s=seq[p-3:p+5]
    elif strand=="-":
        if p<4 or p+4>L:return None
        s=revcomp(seq[p-4:p+4])
    else:return None
    if len(s)!=8 or s[3:5]!="CG" or any(x not in BASES for x in s):return None
    return s
def covband(c):
    if 10<=c<15:return "10-14x"
    if 15<=c<20:return "15-19x"
    if 20<=c<30:return "20-29x"
    if c>=30:return "30plus"
    return None
def ask_file(msg):
    while True:
        p=input(msg).strip().strip('"')
        if os.path.isfile(p):return p
        print("File not found; try again.")

def main():
    print("\nSite-level CpG modeling dataset preparer v1\n")
    genome=ask_file("GRCh38 FASTA (.fa.gz):\n> ");gtf=ask_file("GENCODE GTF (.gtf.gz):\n> ");cgi=ask_file("UCSC cpgIslandExt.txt.gz:\n> ")
    try:k=int(input("Maximum sampled sites per sample x coverage x compartment stratum [5000]: ").strip() or "5000")
    except:k=5000
    entries=[]
    print("Add CpG bedMethyl files. Press Enter when finished.")
    while True:
        p=input("CpG bedMethyl (.bed.gz), or Enter to finish:\n> ").strip().strip('"')
        if not p:break
        if not os.path.isfile(p):print("File not found");continue
        label=input("Sample label: ").strip()
        if not label:label=Path(p).name.split(".")[0]
        entries.append((p,re.sub(r"[^A-Za-z0-9_-]+","_",label)))
    if not entries:return
    out=Path.home()/"Downloads"/"CpG_site_level_model_data";out.mkdir(parents=True,exist_ok=True)
    gidx=load_gtf(gtf);cidx=load_cgi(cgi)
    canonical={f"chr{i}" for i in range(1,23)}|{"chrX","chrY","chrM"}
    cache=Path(tempfile.mkdtemp(prefix="grch38_model_"));found=set();print("Caching genome...")
    for c,seq in fasta_records(genome):
        if c in canonical:
            with gzip.open(cache/f"{c}.seq.gz","wt",compresslevel=1) as f:f.write(seq)
            found.add(c)
    rng=random.Random(20260815)
    reservoirs=defaultdict(list); totals=defaultdict(int); qc=defaultdict(int)
    try:
        for path,label in entries:
            print("Scanning",label)
            current=None;seq=None
            for chrom,p,strand,cov,meth in bed_rows(path):
                qc[(label,"total")]+=1
                if chrom!=current:
                    current=chrom;cp=cache/f"{chrom}.seq.gz"
                    if cp.exists():
                        with gzip.open(cp,"rt") as f:seq=f.read().strip().upper()
                    else:seq=None
                    print(" ",chrom)
                if seq is None:qc[(label,"missing_chrom")]+=1;continue
                band=covband(cov)
                if band is None:qc[(label,"below10")]+=1;continue
                ctx=cpg_context(seq,p,strand)
                if ctx is None:qc[(label,"bad_context")]+=1;continue
                comp=primary(gidx,cidx,chrom,p)
                key=(label,band,comp);totals[key]+=1
                mreads=int(round(cov*meth/100.0));mreads=max(0,min(cov,mreads));ureads=cov-mreads
                row=(label,chrom,p,strand,cov,meth,mreads,ureads,band,comp,ctx[0],ctx[1],ctx[2],ctx[5],ctx[6],ctx[7],ctx[1:7],f"{chrom}:{p//1000000}")
                r=reservoirs[key]
                n=totals[key]
                if len(r)<k:r.append(row)
                else:
                    j=rng.randrange(n)
                    if j<k:r[j]=row
                qc[(label,"kept")]+=1
    finally:shutil.rmtree(cache,ignore_errors=True)
    sampled={key:len(v) for key,v in reservoirs.items()}
    fn=out/"cpg_site_level_sample.tsv.gz"
    header=["sample","chrom","start","strand","coverage","methylation_percent","methylated_reads","unmethylated_reads","coverage_band","compartment","m3","m2","m1","p1","p2","p3","sixmer","block_1Mb","stratum_total","stratum_sampled","sampling_weight"]
    with gzip.open(fn,"wt",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(header)
        for key in sorted(reservoirs):
            total=totals[key];sn=sampled[key];wt=total/sn if sn else 1.0
            for row in reservoirs[key]:w.writerow(list(row)+[total,sn,wt])
    with open(out/"stratum_counts.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["sample","coverage_band","compartment","total_sites","sampled_sites","sampling_weight"])
        for key in sorted(totals):
            sn=sampled[key];w.writerow(list(key)+[totals[key],sn,totals[key]/sn if sn else ""])
    with open(out/"QC.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["sample","metric","value"])
        for (s,m),v in sorted(qc.items()):w.writerow([s,m,v])
    (out/"README.txt").write_text("Upload cpg_site_level_sample.tsv.gz plus stratum_counts.csv and QC.csv for site-level multivariable modeling. Sampling is deterministic and balanced within sample x coverage-band x genomic-compartment strata.\n")
    print("\nDone. Upload the folder:",out)
if __name__=="__main__":main()
