Minimal example of a partial [FastAPI](https://fastapi.tiangolo.com/)-generated API for querying a BGZF-compressed & faidx-indexed FASTA file using [pysam](https://pysam.readthedocs.io/).

Example:

```
$ docker-compose up -d --build
...
$ curl 'http://localhost/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/references'
{"references":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl 'http://localhost/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/glyma.Wm82.gnm2.Gm01:1-100'
{"sequence":"GTTTGGTGTTTGGGTTTTAGGTTTTAGGTTTTAGGTTTTACGGTTTAGGGTTTATGGTTTATGGTTTAGGGTTTAGGGTTAGGAAATAATTTGGGTCTT"}
$curl http://localhost/gff/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz/contigs
{"contigs":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl http://localhost/gff/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz/glyma.Wm82.gnm2.Gm01:1-100000
[{"contig":"glyma.Wm82.gnm2.Gm01","feature":"gene","source":"phytozomev10","start":27354,"end":28320,"score":null,"strand":"-","frame":null,"attributes":"ID=glyma.Wm82.gnm2.ann1.Glyma.01G000100;...
```

See also http://localhost/docs for the FastAPI [Interactive API docs](https://fastapi.tiangolo.com/#interactive-api-docs)
