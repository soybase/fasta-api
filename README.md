Minimal example of a partial [FastAPI](https://fastapi.tiangolo.com/)-generated API for querying a BGZF-compressed & faidx-indexed FASTA file using [pysam](https://pysam.readthedocs.io/).

Example:

```
% make install # one-time environment setup
...
% make
...
$ curl 'http://localhost:8000/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/references'
{"references":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl 'http://localhost:8000/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/glyma.Wm82.gnm2.Gm01:1-100'
{"sequence":"GTTTGGTGTTTGGGTTTTAGGTTTTAGGTTTTAGGTTTTACGGTTTAGGGTTTATGGTTTATGGTTTAGGGTTTAGGGTTAGGAAATAATTTGGGTCTT"}
$ curl http://localhost:8000/gff/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz/contigs
{"contigs":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl http://localhost:8000/gff/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz/glyma.Wm82.gnm2.Gm01:1-100000
[{"contig":"glyma.Wm82.gnm2.Gm01","feature":"gene","source":"phytozomev10","start":27354,"end":28320,"score":null,"strand":"-","frame":null,"attributes":"ID=glyma.Wm82.gnm2.ann1.Glyma.01G000100;...
$ curl http://localhost:8000/vcf/glyma.Wm82.gnm1.div.ContrerasSoto_Mora_2017.SNPs.vcf.gz/contigs
{"contigs":["scaffold_148","scaffold_2079","scaffold_639","scaffold_648","scaffold_1961","scaffold_1902","scaffold_1416","scaffold_1649","scaffold_2267",...
$ curl http://localhost:8000/vcf/glyma.Wm82.gnm1.div.ContrerasSoto_Mora_2017.SNPs.vcf.gz/glyma.Wm82.gnm1.Gm16:1-100000
[{"chrom":"glyma.Wm82.gnm1.Gm16","pos":35846,"id":"M4191","ref":"A","alts":["G"],"qual":null,"filter":[],"info":[],"format":["GT"],"samples":["ANTA","A6001-RR"...
```

See also http://localhost:8000/docs for the FastAPI [Interactive API docs](https://fastapi.tiangolo.com/#interactive-api-docs)
