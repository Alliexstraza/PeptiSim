# Usar imagem oficial Python 3.9 (compatível com suas libs)
FROM python:3.9-slim

# Evitar perguntas durante instalação de pacotes
ENV DEBIAN_FRONTEND=noninteractive

# Instalar dependências de sistema necessárias (exemplo comum para libsbml e compilação)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libsbml5 \
    libsbml-dev \
    libxml2-dev \
    libxslt-dev \
    && rm -rf /var/lib/apt/lists/*

# Criar pasta do app
WORKDIR /app

# Copiar requirements e instalar dependências Python
COPY requirements.txt .

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Copiar todo o código
COPY . .

# Expõe a porta que o Streamlit usa
EXPOSE 8501

# Comando para rodar o app
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
