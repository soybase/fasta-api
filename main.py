# https://pysam.readthedocs.io/en/latest/api.html#fasta-files
import pysam
from fastapi import FastAPI

app = FastAPI()

@app.get("/ref/wm82.gnm2/{seqid}:{start}-{end}")
def wm82_gnm2(seqid: str, start: int, end: int):
  with pysam.FastaFile("https://www.soybase.org/data/v2/Glycine/max/genomes/Wm82.gnm2.DTC4/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz") as f:
    seq = f.fetch(reference="glyma.Wm82.gnm2." + seqid, start = start, end = end)
    return { "sequence" : seq }
