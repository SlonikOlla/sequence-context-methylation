#!/usr/bin/env python3
"""
Annotate GSE70090 site-level reservoir samples with hg19 genomic compartments.

Inputs:
  1) GENCODE v19 comprehensive CHR GTF (.gtf.gz)
  2) UCSC hg19 cpgIslandExt.txt.gz
  3) GSE70090_liver_context_results folder created previously

Outputs:
  GSE70090_oxBS_site_level_annotated.csv.gz
  GSE70090_BS_site_level_annotated.csv.gz
  annotation_QC.csv

Each output retains the sampled CpG-level methylated-read and coverage counts,
local sequence, coverage band, paired patient/status information, a 1-Mb
genomic block, and genomic-compartment annotations.

No external Python packages are required.
"""

import csv, gzip, os, re, sys, bisect
from pathlib import Path
from collections import defaultdict, Counter

def clean(s):
    return s.strip().strip('"').strip("'")

def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if str(path).lower().endswith(".gz") else open(path, "rt", encoding="utf-8", errors="ignore")

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    out=[list(intervals[0])]
    for s,e in intervals[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1]=e
        else:
            out.append([s,e])
    return [(s,e) for s,e in out]

def build_index(raw):
    idx={}
    for key,ints in raw.items():
        merged=merge_intervals(ints)
        starts=[x[0] for x in merged]
        idx[key]=(starts,merged)
    return idx

def contains(idx,key,pos0):
    x=idx.get(key)
    if not x: return False
    starts,ints=x
    i=bisect.bisect_right(starts,pos0)-1
    return i>=0 and ints[i][0] <= pos0 < ints[i][1]

def parse_attrs(s):
    d={}
    for item in s.strip().split(";"):
        item=item.strip()
        if not item: continue
        if " " in item:
            k,v=item.split(" ",1)
            d[k]=v.strip().strip('"')
    return d

def load_gtf(path):
    print("\nReading GENCODE v19 GTF...")
    raw=defaultdict(list)
    counts=Counter()
    with open_text(path) as f:
        for i,line in enumerate(f,1):
            if i%1_000_000==0:
                print(f"  {i:,} GTF lines")
            if not line or line.startswith("#"): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<9: continue
            chrom,feature=p[0],p[2]
            if not chrom.startswith("chr"):
                chrom="chr"+chrom
            if chrom=="chrMT": chrom="chrM"
            try:
                start0=int(p[3])-1
                end0=int(p[4])
            except:
                continue
            strand=p[6]
            attrs=parse_attrs(p[8])
            gene_type=attrs.get("gene_type",attrs.get("gene_biotype",""))

            if feature=="gene":
                raw[(chrom,"gene")].append((start0,end0))
                if gene_type=="protein_coding":
                    raw[(chrom,"pc_gene")].append((start0,end0))
                else:
                    raw[(chrom,"nc_gene")].append((start0,end0))
                # Promoter: TSS -2000/+500, strand-aware.
                if strand=="+":
                    tss=start0
                    ps=max(0,tss-2000); pe=tss+500
                else:
                    tss=end0-1
                    ps=max(0,tss-500); pe=tss+2000
                raw[(chrom,"promoter")].append((ps,pe))
            elif feature=="exon":
                raw[(chrom,"exon")].append((start0,end0))
            elif feature=="CDS":
                raw[(chrom,"CDS")].append((start0,end0))
            counts[feature]+=1
    idx=build_index(raw)
    print("GENCODE features indexed:", dict(counts))
    return idx

def load_cgi(path):
    print("\nReading UCSC hg19 CpG islands...")
    raw=defaultdict(list)
    n=0
    with open_text(path) as f:
        for line in f:
            if not line.strip(): continue
            p=line.rstrip("\n").split("\t")
            # UCSC table: bin, chrom, chromStart, chromEnd, ...
            if len(p)<4: continue
            try:
                if p[0].isdigit() and p[1].startswith("chr"):
                    chrom=p[1]; s=int(p[2]); e=int(p[3])
                elif p[0].startswith("chr"):
                    chrom=p[0]; s=int(p[1]); e=int(p[2])
                else:
                    continue
            except:
                continue
            raw[(chrom,"CGI")].append((s,e)); n+=1
    idx=build_index(raw)
    print(f"Indexed {n:,} CpG islands.")
    return idx

def primary(flags):
    if flags["CpG_island"]: return "CpG_island"
    if flags["promoter"]: return "promoter_non_CGI"
    if flags["CDS"]: return "CDS_nonpromoter_nonCGI"
    if flags["exon"]: return "exon_noncoding_nonpromoter_nonCGI"
    if flags["genic"]: return "intron_or_other_genic_nonpromoter_nonCGI"
    return "intergenic_nonCGI"

def main():
    print("\nGSE70090 hg19 site-level compartment annotator v1\n")
    gtf=clean(input("GENCODE v19 comprehensive CHR GTF (.gtf.gz):\n> "))
    cgi=clean(input("\nUCSC hg19 cpgIslandExt.txt.gz:\n> "))
    folder=clean(input("\nGSE70090_liver_context_results folder:\n> "))

    for p in (gtf,cgi,folder):
        if not os.path.exists(p):
            sys.exit(f"Not found: {p}")

    gidx=load_gtf(gtf)
    cidx=load_cgi(cgi)

    folder=Path(folder)
    files=sorted(folder.glob("P*_site_sample.csv"))
    if len(files)!=16:
        print(f"\nWARNING: expected 16 site-sample files, found {len(files)}")
        for f in files: print(" ",f.name)
        if not files: return

    outpaths={
        "oxBS":folder/"GSE70090_oxBS_site_level_annotated.csv.gz",
        "BS":folder/"GSE70090_BS_site_level_annotated.csv.gz"
    }

    outs={}
    writers={}
    header_out=None
    qc=Counter()

    try:
        for assay,path in outpaths.items():
            outs[assay]=gzip.open(path,"wt",newline="",encoding="utf-8")

        for fi,fpath in enumerate(files,1):
            print(f"\n[{fi}/{len(files)}] {fpath.name}")
            with open(fpath,newline="",encoding="utf-8") as f:
                r=csv.DictReader(f)
                if header_out is None:
                    header_out=list(r.fieldnames)+[
                        "position0","block_1Mb",
                        "CpG_island","promoter","genic","exon","CDS",
                        "protein_coding_gene","noncoding_gene",
                        "primary_compartment"
                    ]
                    for assay in ["oxBS","BS"]:
                        writers[assay]=csv.DictWriter(outs[assay],fieldnames=header_out)
                        writers[assay].writeheader()

                n=0
                for row in r:
                    n+=1
                    chrom=row["chrom"]
                    # raw coordinate in GSE70090 files is 1-based CpG C coordinate
                    pos0=int(row["coordinate_raw"])-1
                    flags={
                        "CpG_island":contains(cidx,(chrom,"CGI"),pos0),
                        "promoter":contains(gidx,(chrom,"promoter"),pos0),
                        "genic":contains(gidx,(chrom,"gene"),pos0),
                        "exon":contains(gidx,(chrom,"exon"),pos0),
                        "CDS":contains(gidx,(chrom,"CDS"),pos0),
                        "protein_coding_gene":contains(gidx,(chrom,"pc_gene"),pos0),
                        "noncoding_gene":contains(gidx,(chrom,"nc_gene"),pos0),
                    }
                    row["position0"]=pos0
                    row["block_1Mb"]=f"{chrom}:{pos0//1_000_000}"
                    for k,v in flags.items():
                        row[k]=1 if v else 0
                    row["primary_compartment"]=primary(flags)
                    assay=row["assay"]
                    writers[assay].writerow(row)
                    qc[f"rows_{assay}"]+=1
                    qc[f"primary_{row['primary_compartment']}"]+=1
                    if n%100_000==0:
                        print(f"  {n:,} rows annotated")
    finally:
        for h in outs.values():
            h.close()

    with open(folder/"annotation_QC.csv","w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["metric","value"])
        for k,v in sorted(qc.items()): w.writerow([k,v])

    print("\nDONE.")
    print("Created:")
    print(outpaths["oxBS"])
    print(outpaths["BS"])
    print("\nFor the primary cancer model, ZIP and upload:")
    print(outpaths["oxBS"])
    print(folder/"annotation_QC.csv")
    print("\nBS can be uploaded afterward as confirmation.")

if __name__=="__main__":
    main()
