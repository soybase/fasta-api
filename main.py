# https://pysam.readthedocs.io/en/latest/api.html#fasta-files
import json
import pysam
from fastapi import FastAPI

BASE_URL = "https://www.soybase.org/data/"

app = FastAPI()

with open("data.json") as f:
    data = json.load(f)

@app.get("/fasta/{fasta}/{seqid}:{start}-{end}")
def wm82_gnm2(fasta: str, seqid: str, start: int, end: int):
  with pysam.FastaFile(BASE_URL + data[fasta] + fasta) as f:
    seq = f.fetch(reference=seqid, start = start, end = end)
    return { "sequence" : seq }

@app.get("/fasta/{fasta}/references")
def wm82_gnm2(fasta: str):
  with pysam.FastaFile(BASE_URL + data[fasta] + fasta) as f:
    return { "references": f.references }
