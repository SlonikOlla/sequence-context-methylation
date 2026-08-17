#!/usr/bin/env python3
"""
Annotate ENCODE CpG bedMethyl with genomic compartments and sequence context.

Inputs:
  - GRCh38 primary-assembly FASTA (.fa.gz or .fa)
  - ENCODE CpG bedMethyl (.bed.gz)
  - GENCODE v50 comprehensive GTF, chromosome version (.gtf.gz)
  - UCSC hg38 cpgIslandExt.txt.gz

Outputs compact summaries only:
  - QC
  - compartment overall methylation
  - compartment x position/base
  - compartment x methylation-bin x position/base
  - compartment x 6-mer
  - compartment x 8-mer
  - selected motif summary
  - chromosome summary

Primary mutually exclusive compartments:
  CpG_island
  promoter_non_CGI
  CDS_nonpromoter_nonCGI
  exon_noncoding_nonpromoter_nonCGI
  intron_nonpromoter_nonCGI
  other_genic_nonpromoter_nonCGI
  intergenic_nonCGI

Also outputs binary flags for:
  genic/intergenic, protein_coding_gene/noncoding_gene, promoter, exon,
  CDS, CpG_island.

Promoter definition: TSS -2000 bp to +500 bp, strand aware.
Coverage threshold default: 10.
"""
import csv, gzip, os, re, shutil, tempfile, bisect
from pathlib import Path
from collections import defaultdict

BASES=set("ACGT")
COMP=str.maketrans("ACGTNacgtn","TGCANtgcan")

def revcomp(s): return s.translate(COMP)[::-1].upper()
def op_text(p,mode="rt"):
    return gzip.open(p,mode) if str(p).lower().endswith(".gz") else open(p,mode)

def fasta_records(path):
    with op_text(path) as f:
        name=None; seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None: yield name,"".join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line.strip())
        if name is not None: yield name,"".join(seq).upper()

def merge_intervals(iv):
    if not iv: return []
    iv=sorted(iv)
    out=[list(iv[0])]
    for s,e in iv[1:]:
        if s<=out[-1][1]:
            if e>out[-1][1]: out[-1][1]=e
        else: out.append([s,e])
    return [(s,e) for s,e in out]

def point_in(merged, starts, p):
    if not merged: return False
    i=bisect.bisect_right(starts,p)-1
    return i>=0 and p<merged[i][1]

def parse_attrs(s):
    d={}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"',s):
        d[m.group(1)]=m.group(2)
    return d

def is_noncoding_gene_type(gt):
    if not gt: return False
    return gt not in {"protein_coding"}

def load_gtf(gtf):
    raw=defaultdict(lambda: defaultdict(list))
    print("Reading GENCODE GTF...")
    with op_text(gtf) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<9: continue
            chrom,_,feature,start,end,_,strand,_,attrs=p
            if not chrom.startswith("chr"): continue
            s=int(start)-1; e=int(end)   # GTF -> 0-based half-open
            a=parse_attrs(attrs)
            gt=a.get("gene_type") or a.get("gene_biotype") or ""
            if feature=="gene":
                raw[chrom]["gene"].append((s,e))
                if gt=="protein_coding": raw[chrom]["pc_gene"].append((s,e))
                else: raw[chrom]["nc_gene"].append((s,e))
                # strand-aware promoter around TSS
                tss=s if strand=="+" else e-1
                ps=max(0,tss-2000) if strand=="+" else max(0,tss-500)
                pe=tss+501 if strand=="+" else tss+2001
                raw[chrom]["promoter"].append((ps,pe))
            elif feature=="exon":
                raw[chrom]["exon"].append((s,e))
            elif feature=="CDS":
                raw[chrom]["CDS"].append((s,e))
    idx={}
    for chrom,types in raw.items():
        idx[chrom]={}
        for typ,iv in types.items():
            m=merge_intervals(iv)
            idx[chrom][typ]=(m,[x[0] for x in m])
        print(" ",chrom, {k:len(v[0]) for k,v in idx[chrom].items()})
    return idx

def load_cgi(path):
    raw=defaultdict(list)
    print("Reading UCSC CpG islands...")
    with op_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            p=line.rstrip("\n").split("\t")
            # UCSC cpgIslandExt.txt.gz begins: bin, chrom, chromStart, chromEnd...
            try:
                if p[0].isdigit() and len(p)>=4 and p[1].startswith("chr"):
                    chrom=p[1]; s=int(p[2]); e=int(p[3])
                elif p[0].startswith("chr") and len(p)>=3:
                    chrom=p[0]; s=int(p[1]); e=int(p[2])
                else: continue
            except: continue
            raw[chrom].append((s,e))
    idx={}
    for chrom,iv in raw.items():
        m=merge_intervals(iv)
        idx[chrom]=(m,[x[0] for x in m])
    print(" CpG islands loaded:",sum(len(v[0]) for v in idx.values()))
    return idx

def getflag(idx,chrom,typ,p):
    x=idx.get(chrom,{}).get(typ)
    return point_in(x[0],x[1],p) if x else False

def getcgi(idx,chrom,p):
    x=idx.get(chrom)
    return point_in(x[0],x[1],p) if x else False

def primary_compartment(flags):
    if flags["CpG_island"]: return "CpG_island"
    if flags["promoter"]: return "promoter_non_CGI"
    if flags["CDS"]: return "CDS_nonpromoter_nonCGI"
    if flags["exon"]: return "exon_noncoding_nonpromoter_nonCGI"
    if flags["genic"]: return "intron_or_other_genic_nonpromoter_nonCGI"
    return "intergenic_nonCGI"

def bed_rows(path):
    with op_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<11: p=line.rstrip("\n").split()
            if len(p)<11: continue
            try:
                yield p[0],int(p[1]),p[5],float(p[9]),float(p[10])
            except: continue

def build_cache(genome,needed,cache):
    found=set()
    print("Building chromosome sequence cache...")
    for chrom,seq in fasta_records(genome):
        if chrom in needed:
            with gzip.open(cache/f"{chrom}.seq.gz","wt",compresslevel=1) as f:f.write(seq)
            found.add(chrom); print(" ",chrom,len(seq))
    return found

def loadseq(cache,chrom):
    p=cache/f"{chrom}.seq.gz"
    if not p.exists(): return None
    with gzip.open(p,"rt") as f:return f.read().strip().upper()

def cpg_context(seq,p,strand):
    L=len(seq)
    if strand=="+":
        if p<3 or p+5>L:return None
        s=seq[p-3:p+5]
    elif strand=="-":
        if p<4 or p+4>L:return None
        s=revcomp(seq[p-4:p+4])
    else:
        if p>=3 and p+5<=L and seq[p:p+2]=="CG": s=seq[p-3:p+5]
        elif p>=4 and p+4<=L and seq[p-1:p+1]=="CG": s=revcomp(seq[p-4:p+4])
        else:return None
    if len(s)!=8 or s[3:5]!="CG" or any(c not in BASES for c in s):return None
    return s

class Stat:
    __slots__=("n","s","ss","cov")
    def __init__(self):self.n=0;self.s=0.;self.ss=0.;self.cov=0.
    def add(self,m,c):self.n+=1;self.s+=m;self.ss+=m*m;self.cov+=c
    def row(self):
        mean=self.s/self.n
        var=max(0.,(self.ss-self.s*self.s/self.n)/(self.n-1)) if self.n>1 else 0.
        return [self.n,mean,var**0.5,self.cov/self.n]

def main():
    print("\nCpG genomic-compartment reducer v2\n")
    def ask(msg):
        while True:
            p=input(msg).strip().strip('"')
            if os.path.isfile(p):return p
            print("File not found; try again.")
    genome=ask("GRCh38 primary-assembly FASTA (.fa.gz):\n> ")
    bed=ask("ENCODE CpG bedMethyl (.bed.gz):\n> ")
    gtf=ask("GENCODE v50 comprehensive CHR GTF (.gtf.gz):\n> ")
    cgi=ask("UCSC hg38 cpgIslandExt.txt.gz:\n> ")
    try:mincov=int(input("Minimum coverage [10]: ").strip() or "10")
    except:mincov=10

    sample=input("Short sample label (e.g. Liver53F or HepG2): ").strip() or Path(bed).name.split(".")[0]
    sample=re.sub(r"[^A-Za-z0-9_-]+","_",sample)
    out=Path.home()/"Downloads"/f"WGBS_compartment_results_{sample}"
    out.mkdir(parents=True,exist_ok=True)
    test=out/"_write_test.tmp";test.write_text("ok");test.unlink()
    print("Output:",out)

    gidx=load_gtf(gtf)
    cidx=load_cgi(cgi)

    # collect BED chromosomes
    needed=set()
    print("Reading BED chromosome names...")
    for chrom,_,_,_,_ in bed_rows(bed):needed.add(chrom)

    cache=Path(tempfile.mkdtemp(prefix="grch38_ctx_"))
    found=build_cache(genome,needed,cache)

    overall=defaultdict(Stat)                # category
    posstat=defaultdict(Stat)                # category,pos,base
    binpos=defaultdict(int)                  # category,bin,pos,base
    sixstat=defaultdict(Stat)                # category,sixmer
    eightstat=defaultdict(Stat)              # primary,8mer
    chromstat=defaultdict(Stat)
    qc=defaultdict(int)
    current=None; seq=None
    positions=[("-3",0),("-2",1),("-1",2),("+1",5),("+2",6),("+3",7)]

    print("\nScanning WGBS and annotating compartments...")
    try:
        for chrom,p,strand,cov,meth in bed_rows(bed):
            qc["total_rows"]+=1
            if chrom!=current:
                current=chrom;seq=loadseq(cache,chrom) if chrom in found else None
                print(" ",chrom)
            if seq is None:
                qc["missing_chrom"]+=1;continue
            if cov<mincov:
                qc["below_coverage"]+=1;continue
            ctx=cpg_context(seq,p,strand)
            if ctx is None:
                qc["bad_context"]+=1;continue
            qc["kept"]+=1

            flags={
                "CpG_island":getcgi(cidx,chrom,p),
                "promoter":getflag(gidx,chrom,"promoter",p),
                "genic":getflag(gidx,chrom,"gene",p),
                "exon":getflag(gidx,chrom,"exon",p),
                "CDS":getflag(gidx,chrom,"CDS",p),
                "protein_coding_gene":getflag(gidx,chrom,"pc_gene",p),
                "noncoding_gene":getflag(gidx,chrom,"nc_gene",p),
            }
            flags["intergenic"]=not flags["genic"]
            primary=primary_compartment(flags)

            categories=["ALL","primary:"+primary]
            for k,v in flags.items():
                if v: categories.append("flag:"+k)

            six=ctx[1:7]
            b=min(9,int(meth//10))
            binlab=f"{b*10}-{100 if b==9 else (b+1)*10}%"
            for cat in categories:
                overall[cat].add(meth,cov)
                sixstat[(cat,six)].add(meth,cov)
                for pn,ix in positions:
                    base=ctx[ix]
                    posstat[(cat,pn,base)].add(meth,cov)
                    binpos[(cat,binlab,pn,base)]+=1
            eightstat[(primary,ctx)].add(meth,cov)
            chromstat[chrom].add(meth,cov)
    finally:
        shutil.rmtree(cache,ignore_errors=True)

    def write_stats(filename,header,items):
        with open(out/filename,"w",newline="") as f:
            w=csv.writer(f);w.writerow(header)
            for keys,st in items:
                if not isinstance(keys,tuple):keys=(keys,)
                w.writerow(list(keys)+st.row())

    with open(out/"QC.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["metric","value"])
        w.writerow(["coverage_threshold",mincov])
        for k in ["total_rows","kept","below_coverage","bad_context","missing_chrom"]:w.writerow([k,qc[k]])

    write_stats("compartment_overall.csv",
        ["category","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(overall.items()))
    write_stats("compartment_position.csv",
        ["category","position","base","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(posstat.items()))
    with open(out/"compartment_bin_position.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["category","methylation_bin","position","base","site_count"])
        for k,n in sorted(binpos.items()):w.writerow(list(k)+[n])
    write_stats("compartment_6mer.csv",
        ["category","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(sixstat.items()))
    write_stats("primary_compartment_8mer.csv",
        ["primary_compartment","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(eightstat.items()))
    write_stats("chromosome_summary.csv",
        ["chromosome","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"],
        sorted(chromstat.items()))

    # selected motifs
    selected={"TACGTA","TCCGTA","TACGGA","GGCGCG","TGCGGT","AGCGTA","AACGTT"}
    with open(out/"selected_motifs_by_compartment.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["category","motif","site_count","mean_methylation_percent","sd_methylation_percent","mean_coverage"])
        for (cat,m),st in sorted(sixstat.items()):
            if m in selected:w.writerow([cat,m]+st.row())

    with open(out/"README.txt","w") as f:
        f.write("CpG genomic compartment analysis\n")
        f.write("Promoter: strand-aware TSS -2000/+500 bp.\n")
        f.write("Primary hierarchy: CpG island > promoter > CDS > exon > genic > intergenic.\n")
        f.write("Binary flags are also summarized independently.\n")
        f.write(f"Coverage threshold: {mincov}x\n")
        f.write("ZIP this folder and upload it to ChatGPT.\n")

    print("\nDONE")
    print("Kept CpGs:",f"{qc['kept']:,}")
    print("Results:",out)
    print("ZIP WGBS_compartment_results and upload it.")

if __name__=="__main__":
    main()
