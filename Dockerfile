# Usamos uma imagem Python oficial, slim para ser leve
FROM python:3.9-slim

# Atualiza o apt e instala a libsbml (dependência nativa do python-libsbml)
RUN apt-get update && apt-get install -y libsbml5 && rm -rf /var/lib/apt/lists/*

# Define o diretório de trabalho dentro do container
WORKDIR /app

# Copia todos os arquivos do seu projeto para o container
COPY . /app

# Atualiza pip e instala as dependências Python
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Expõe a porta 8501 (default do Streamlit)
EXPOSE 8501

# Comando para rodar o app Streamlit, ouvindo em todas interfaces e porta 8501
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
