# REFERENCES
#     https://fastapi.tiangolo.com/deployment/docker/
FROM ubuntu:20.04

RUN apt update && apt install -y --no-install-recommends \
  python3-pysam \
  python3-pip \
  && pip3 install --no-cache-dir fastapi==0.65.1 uvicorn==0.13.4 \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY data.json .
COPY main.py .

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]
