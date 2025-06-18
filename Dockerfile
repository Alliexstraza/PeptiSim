FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake libsbml5 libsbml-dev libxml2-dev libxslt-dev curl \
    && rm -rf /var/lib/apt/lists/*

# Instala o solver HiGHS — versão 1.9.0
RUN HIGHS_VER=1.9.0 && \
    curl -L -o highs.tar.gz \
      https://github.com/ERGO-Code/HiGHS/releases/download/v${HIGHS_VER}/highs-x86_64-Linux.tar.gz && \
    file highs.tar.gz && \
    tar -xzf highs.tar.gz && \
    mv highs/bin/highs /usr/local/bin/ && \
    rm -rf highs highs.tar.gz

WORKDIR /app

COPY requirements.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
