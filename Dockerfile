FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libsbml5 \
    libsbml-dev \
    libxml2-dev \
    libxslt-dev \
    curl \
    tar \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Baixa o arquivo
RUN curl -L -o highs.tar.gz https://github.com/Alliexstraza/PeptiSim/raw/refs/heads/main/binaries/HiGHS-1.9.0.tar.gz

# Lista para garantir que baixou
RUN ls -lh highs.tar.gz

# Extrai o arquivo
RUN tar -xzf highs.tar.gz

# Lista a pasta extraída
RUN ls -l HiGHS-1.9.0/
RUN ls -l HiGHS-1.9.0/bin/

# Move o binário para /usr/local/bin
RUN mv HiGHS-1.9.0/bin/highs /usr/local/bin/

# Remove arquivos temporários
RUN rm -rf HiGHS-1.9.0 highs.tar.gz

# Continua com o resto do seu Dockerfile...

WORKDIR /app

COPY requirements.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
