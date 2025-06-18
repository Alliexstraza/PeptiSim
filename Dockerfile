# Usa Python 3.9 slim como base (leve e compatível com COBRApy)
FROM python:3.9-slim

# Impede prompts interativos durante a instalação
ENV DEBIAN_FRONTEND=noninteractive

# Instala dependências de sistema (incluindo o básico pro solver HiGHS)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libsbml5 \
    libsbml-dev \
    libxml2-dev \
    libxslt-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Instala o solver HiGHS (baixado do seu repositório no GitHub)
RUN curl -L -o highs.tar.gz \
    https://github.com/Alliexstraza/PeptiSim/raw/refs/heads/main/binaries/HiGHS-1.9.0.tar.gz && \
    tar -xzf highs.tar.gz && \
    mv HiGHS-1.9.0/bin/highs /usr/local/bin/ && \
    rm -rf HiGHS-1.9.0 highs.tar.gz

# Define diretório de trabalho
WORKDIR /app

# Instala dependências do Python
COPY requirements.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt

# Copia os arquivos do app
COPY . .

# Expõe a porta padrão do Streamlit
EXPOSE 8501

# Comando padrão para rodar o app no Render
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
