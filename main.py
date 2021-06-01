# https://pysam.readthedocs.io/en/latest/api.html#fasta-files
import json
import pysam
from fastapi import FastAPI

BASE_URL = "https://www.soybase.org/data/"

app = FastAPI()

with open("data.json") as f:
    data = json.load(f)

@app.get("/fasta/{fasta}/{seqid}:{start}-{end}")
def fasta_range(fasta: str, seqid: str, start: int, end: int):
    seq = pysam.FastaFile(BASE_URL + data[fasta] + fasta).fetch(reference=seqid, start = start, end = end)
    return { "sequence" : seq }

@app.get("/fasta/{fasta}/references")
def fasta_references(fasta: str):
    return { "references": pysam.FastaFile(BASE_URL + data[fasta] + fasta).references }

@app.get("/gff/{gff}/contigs")
def gff_references(gff: str):
    return { "contigs": pysam.TabixFile(BASE_URL + data[gff] + gff).contigs }

@app.get("/gff/{gff}/{seqid}:{start}-{end}")
def gff_features(gff: str, seqid: str, start: int, end: int):
  return [ {"contig": feature.contig,
            "feature": feature.feature,
            "source": feature.source,
            "start": feature.start,
            "end": feature.end,
            "score": feature.score,
            "strand": feature.strand,
            "frame": feature.frame,
            "attributes": list(feature.attributes)} 
            for feature 
            in pysam.TabixFile(BASE_URL + data[gff] + gff).fetch(seqid, start, end, parser=pysam.asGFF3()) ]

@app.get("/vcf/{vcf}/contigs")
def vcf_contigs(vcf: str):
  return { "contigs": list(pysam.VariantFile(BASE_URL + data[vcf] + vcf).header.contigs) }

@app.get("/vcf/{vcf}/{seqid}:{start}-{end}")
def vcf_features(vcf: str, seqid: str, start: int, end: int):
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
            in pysam.VariantFile(BASE_URL + data[vcf] + vcf).fetch(seqid, start, end) ]
