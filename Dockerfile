# Usa Python 3.9 slim (leve e compatível com COBRApy)
FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

# Instala dependências de sistema e o solver HiGHS
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libsbml5 \
    libsbml-dev \
    libxml2-dev \
    libxslt-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Instala o solver HiGHS
RUN curl -L https://github.com/ERGO-Code/HiGHS/releases/download/v1.6.0/highs-linux64.tar.gz | tar xz \
    && mv highs/bin/highs /usr/local/bin/ \
    && rm -r highs

WORKDIR /app

COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
