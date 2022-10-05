# https://pysam.readthedocs.io/en/latest/api.html#fasta-files
import json
import pysam
from fastapi import FastAPI
import urllib
import itertools


app = FastAPI()

with open("data.json") as f:
    data = json.load(f)

@app.get("/fasta/fetch/{seqid}:{start}-{end}/{url:path}")
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

@app.get("/gff/fetch/{seqid}:{start}-{end}/{url:path}")
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

#As itemRgb, blockSizes, and blockStarts columns are rare, their types have not been determined and may change.
@app.get("/bed/fetch/{seqid}:{start}-{end}/{url:path}")
def bed_features(url: str, seqid: str, start: int, end: int):
  bedcols = ('contig', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRGB', 'blockCount', 'blockSizes', 'blockStarts')
  return [{
            "data": dict(itertools.starmap(lambda k,v: (k, int(v) if k in ['start', 'end', 'thickStart', 'thickEnd'] else float(v) if k == 'score' else v), zip(bedcols, feature)))
            } 
            for feature 
            in pysam.TabixFile(urllib.parse.unquote(url)).fetch(seqid, start, end, parser=pysam.asBed()) ]

@app.get("/vcf/contigs/{url:path}")
def vcf_contigs(url: str):
  return { "contigs": list(pysam.VariantFile(urllib.parse.unquote(url)).header.contigs) }

@app.get("/vcf/fetch/{seqid}:{start}-{end}/{url:path}")
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

@app.get("/alignment/references/{url:path}")
def alignment_references(url: str):
    return { "references": pysam.AlignmentFile(urllib.parse.unquote(url)).references }


@app.get("/alignment/unmapped/{url:path}")
def alignment_unmapped(url: str):
    return { "unmapped": pysam.AlignmentFile(urllib.parse.unquote(url)).unmapped }

@app.get("/alignment/nreferences/{url:path}")
def alignment_nreferences(url: str):
    return { "nreferences": pysam.AlignmentFile(urllib.parse.unquote(url)).nreferences }

@app.get("/alignment/nocoordinate/{url:path}")
def alignment_nocoordinate(url: str):
    return { "nocoordinate": pysam.AlignmentFile(urllib.parse.unquote(url)).nocoordinate }

@app.get("/alignment/mapped/{url:path}")
def alignment_mapped(url: str):
    return { "mapped": pysam.AlignmentFile(urllib.parse.unquote(url)).mapped }

@app.get("/alignment/lengths/{url:path}")
def alignment_lengths(url: str):
    return { "lengths": pysam.AlignmentFile(urllib.parse.unquote(url)).lengths }

@app.get("/alignment/index_statistics/{url:path}")
def alignment_index_statistics(url: str):
    return { "index_statistics": pysam.AlignmentFile(urllib.parse.unquote(url)).get_index_statistics() }




@app.get("/alignment/count/{contig}:{start}-{stop}/{url:path}")
def alignment_count(url: str, contig: str, start: int, stop: int):
    count = pysam.AlignmentFile(urllib.parse.unquote(url)).count(contig, start, stop)
    return { "count" : count }

@app.get("/alignment/count_coverage/{contig}:{start}-{stop}/{url:path}")
def alignment_count_coverage(url: str, contig: str, start: int, stop: int):
    count_coverage = pysam.AlignmentFile(urllib.parse.unquote(url)).count_coverage(contig, start, stop)
    return [{"A" : json.dumps([x for x in count_coverage[0]]),
             "B" : json.dumps([x for x in count_coverage[1]]),
             "C" : json.dumps([x for x in count_coverage[2]]),
             "D" : json.dumps([x for x in count_coverage[3]])
                }]


@app.get("/alignment/fetch/{contig}:{start}-{stop}/{url:path}")
def alignment_fetch(url: str, contig: str, start: int, stop: int):
     return [ {"data": feature.to_dict()}
            for feature
            in pysam.AlignmentFile(urllib.parse.unquote(url)).fetch(contig=contig, start = start, stop = stop) ]


@app.get("/alignment/length/{reference}/{url:path}")
def alignment_lengths(reference: str , url: str):
    return { "length": pysam.AlignmentFile(urllib.parse.unquote(url)).get_reference_length(reference) }


    
