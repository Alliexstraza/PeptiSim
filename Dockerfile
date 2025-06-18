FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    tar \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

RUN curl -L -o highs.tar.gz https://github.com/Alliexstraza/PeptiSim/raw/refs/heads/main/binaries/HiGHS-1.9.0.tar.gz

# Lista o conteúdo do arquivo compactado sem extrair
RUN tar -tzf highs.tar.gz
