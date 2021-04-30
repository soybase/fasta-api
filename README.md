Minimal example of a partial [FastAPI](https://fastapi.tiangolo.com/)-generated API for querying a BGZF-compressed & faidx-indexed FASTA file using [pysam](https://pysam.readthedocs.io/).

Example:

```
$ docker-compose up -d --build
...
$ curl 'http://localhost/ref/wm82.gnm2/Gm01:1-100'
{"sequence":"GTTTGGTGTTTGGGTTTTAGGTTTTAGGTTTTAGGTTTTACGGTTTAGGGTTTATGGTTTATGGTTTAGGGTTTAGGGTTAGGAAATAATTTGGGTCTT"}
```

See also http://localhost/docs for the FastAPI [Interactive API docs](https://fastapi.tiangolo.com/#interactive-api-docs)
