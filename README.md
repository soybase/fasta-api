Minimal example of a partial [FastAPI](https://fastapi.tiangolo.com/)-generated API for querying a BGZF-compressed & faidx-indexed FASTA file using [pysam](https://pysam.readthedocs.io/).

Example:

```
$ docker-compose up -d --build
...
$ curl 'http://localhost/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/glyma.Wm82.gnm2.Gm01:1-100'
{"sequence":"GTTTGGTGTTTGGGTTTTAGGTTTTAGGTTTTAGGTTTTACGGTTTAGGGTTTATGGTTTATGGTTTAGGGTTTAGGGTTAGGAAATAATTTGGGTCTT"}
$ curl 'http://localhost/fasta/glyma.Wm82.gnm2.DTC4.genome_main.fna.gz/references'
{"references":["glyma.Wm82.gnm2.Gm01","glyma.Wm82.gnm2.Gm02",...

```

See also http://localhost/docs for the FastAPI [Interactive API docs](https://fastapi.tiangolo.com/#interactive-api-docs)

data.yml

```
cd /usr/local/www/data
find v2 -name '*.gz' | awk -F '/' '{print $NF ":", "/" $0}' > data.yml
```
