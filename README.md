Minimal example of a partial [FastAPI](https://fastapi.tiangolo.com/)-generated API for querying a BGZF-compressed & faidx-indexed FASTA file using [pysam](https://pysam.readthedocs.io/).

Adapted to run as an Azure function app from the example at https://github.com/Azure-Samples/fastapi-on-azure-functions/ , omitting the Azure Bicep IaC components (instead assuming Azure function app will be created via web UI or Azure CLI)

## Local Development

```
make install
make # or "make test"
```

## Deploy

Create an Azure function app (e.g., in the Azure web portal) called *fasta-api*, then:

```
make login
make publish
make publish # NOTE: for some reason, first deployment doesn't complete, and this needs to be run twice???
```

## Examples

Assuming application endpoint is https://fasta-api.azurewebsites.net

```
% make install # one-time environment setup
...
% make
...
$ curl https://fasta-api.azurewebsites.net/fasta/references/https://www.soybase.org/data/v2/Glycine/max/genomes/Wm82.gnm2.DTC4/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz
{"references":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl https://fasta-api.azurewebsites.net/fasta/fetch/glyma.Wm82.gnm2.Gm01:1-100/https://www.soybase.org/data/v2/Glycine/max/genomes/Wm82.gnm2.DTC4/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz 
{"sequence":"GTTTGGTGTTTGGGTTTTAGGTTTTAGGTTTTAGGTTTTACGGTTTAGGGTTTATGGTTTATGGTTTAGGGTTTAGGGTTAGGAAATAATTTGGGTCTT"}
$ curl https://fasta-api.azurewebsites.net/gff/contigs/https://www.soybase.org/data/v2/Glycine/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz 
{"contigs":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...
$ curl https://fasta-api.azurewebsites.net/gff/fetch/glyma.Wm82.gnm2.Gm01:1-100000/https://www.soybase.org/data/v2/Glycine/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz
[{"contig":"glyma.Wm82.gnm2.Gm01","feature":"gene","source":"phytozomev10","start":27354,"end":28320,"score":null,"strand":"-","frame":null,"attributes":"ID=glyma.Wm82.gnm2.ann1.Glyma.01G000100;...
$ curl https://fasta-api.azurewebsites.net/vcf/contigs/https://www.soybase.org/data/v2/Glycine/max/diversity/Wm82.gnm1.div.ContrerasSoto_Mora_2017/glyma.Wm82.gnm1.div.ContrerasSoto_Mora_2017.SNPs.vcf.gz
{"contigs":["scaffold_148","scaffold_2079","scaffold_639","scaffold_648","scaffold_1961","scaffold_1902","scaffold_1416","scaffold_1649","scaffold_2267",...
$ curl https://fasta-api.azurewebsites.net/vcf/fetch/glyma.Wm82.gnm1.Gm16:1-100000/https://www.soybase.org/data/v2/Glycine/max/diversity/Wm82.gnm1.div.ContrerasSoto_Mora_2017/glyma.Wm82.gnm1.div.ContrerasSoto_Mora_2017.SNPs.vcf.gz
[{"chrom":"glyma.Wm82.gnm1.Gm16","pos":35846,"id":"M4191","ref":"A","alts":["G"],"qual":null,"filter":[],"info":[],"format":["GT"],"samples":["ANTA","A6001-RR"...
$ curl https://fasta-api.azurewebsites.net/bed/fetch/glyma.Wm82_ISU01.gnm2.Gm01:1-100000/https://data.legumeinfo.org/Glycine/max/annotations/Wm82_ISU01.gnm2.ann1.FGFB/glyma.Wm82_ISU01.gnm2.ann1.FGFB.gene_models_main.bed.gz
[{"contig":"glyma.Wm82_ISU01.gnm2.Gm01","start":78502,"end":103594,"name":"glyma.Wm82_ISU01.gnm2.ann1.GmISU01.01G000050.1","score":0.0,"strand":"-"}]
$ curl https://fasta-api.azurewebsites.net/alignment/fetch/aradu.V14167.gnm2.chr01:1-100000/https://data.legumeinfo.org/Arachis/duranensis/genome_alignments/V14167.gnm2.wga.96TT/aradu.V14167.gnm2.x.araca.K10017.gnm1.96TT.bam
[{"name":"araca.K10017.gnm1.chr01","flag":"2048","ref_name":"aradu.V14167.gnm2.chr01","ref_pos":"45008","map_quality":"60","cigar":"191657H159M11I364M10D109M2D685M10I111...
'''



See also https://fasta-api.azurewebsites.net/docs for the FastAPI [Interactive API docs](https://fastapi.tiangolo.com/#interactive-api-docs)
