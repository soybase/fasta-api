# https://pysam.readthedocs.io/en/latest/api.html#fasta-files
import json
import pysam
from fastapi import FastAPI
import urllib


app = FastAPI()

with open("data.json") as f:
    data = json.load(f)

@app.get("/fasta/{seqid}:{start}-{end}/{url:path}")
def fasta_range(url: str, seqid: str, start: int, end: int):
    seq = pysam.FastaFile(urllib.parse.unquote(url)).fetch(reference=seqid, start = start, end = end)
    return { "sequence" : seq }

@app.get("/fasta/references/{url:path}")
def fasta_references(url: str):
    return { "references": pysam.FastaFile(urllib.parse.unquote(url)).references }

@app.get("/fasta/lengths/{url:path}")
def fasta_lengths(url: str):
    return { "lengths": pysam.FastaFile(urllib.parse.unquote(url)).lengths }

@app.get("/fasta/nreferences/{url:path}")
def fasta_nreferences(url: str):
    return { "nreferences": pysam.FastaFile(urllib.parse.unquote(url)).nreferences }

@app.get("/gff/contigs/{url:path}")
def gff_references(url: str):
    return { "contigs": pysam.TabixFile(urllib.parse.unquote(url)).contigs }

@app.get("/gff/{seqid}:{start}-{end}/{url:path}")
def gff_features(url: str, seqid: str, start: int, end: int):
  return [ {"contig": feature.contig,
            "feature": feature.feature,
            "source": feature.source,
            "start": feature.start,
            "end": feature.end,
            "score": feature.score,
            "strand": feature.strand,
            "frame": feature.frame,
            "attributes": dict(a.split("=") for a in feature.attributes.split(";") if a != "")} 
            for feature 
            in pysam.TabixFile(urllib.parse.unquote(url)).fetch(seqid, start, end, parser=pysam.asGFF3()) ]

@app.get("/vcf/contigs/{url:path}")
def vcf_contigs(url: str):
  return { "contigs": list(pysam.VariantFile(urllib.parse.unquote(url)).header.contigs) }

@app.get("/vcf/{seqid}:{start}-{end}/{url:path}")
def vcf_features(url: str, seqid: str, start: int, end: int):
  return [ {"chrom":   feature.chrom,
            "pos":     feature.pos,
            "id":      feature.id,
            "ref":     feature.ref,
            "alts":    feature.alts,
            "qual":    feature.qual,
            "filter":  list(feature.filter),
            "info":    list(feature.info),
            "format":  list(feature.format),
            "samples": list(feature.samples),
            "alleles": feature.alleles}
            for feature 
            in pysam.VariantFile(urllib.parse.unquote(url)).fetch(seqid, start, end) ]
