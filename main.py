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
  with pysam.FastaFile(BASE_URL + data[fasta] + fasta) as f:
    seq = f.fetch(reference=seqid, start = start, end = end)
    return { "sequence" : seq }

@app.get("/fasta/{fasta}/references")
def fasta_references(fasta: str):
  with pysam.FastaFile(BASE_URL + data[fasta] + fasta) as f:
    return { "references": f.references }

@app.get("/gff/{gff}/contigs")
def gff_references(gff: str):
  with pysam.TabixFile(BASE_URL + data[gff] + gff) as f:
    return { "contigs": f.contigs }

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
            "attributes": feature.attributes} 
            for feature 
            in pysam.TabixFile(BASE_URL + data[gff] + gff).fetch(seqid, start, end, parser=pysam.asGFF3()) ]
