FROM python:3.12-slim

WORKDIR /app

# ---------- Variáveis de ambiente ----------
ENV DEBIAN_FRONTEND=noninteractive
ENV GOLANG_VERSION=1.21.5
ENV PATH="/usr/local/go/bin:$PATH"


# ---------- Instala dependências de sistema ----------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    gnupg \
    ca-certificates \
    apt-transport-https \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

# ---------- Instala o Go ----------
RUN wget https://golang.org/dl/go${GOLANG_VERSION}.linux-amd64.tar.gz && \
    tar -C /usr/local -xzf go${GOLANG_VERSION}.linux-amd64.tar.gz && \
    rm go${GOLANG_VERSION}.linux-amd64.tar.gz

# ---------- Instala o Minimap2 ----------
RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && make && cp minimap2 /usr/local/bin && cd .. && rm -rf minimap2

# ---------- Instala o gofasta ----------
RUN go install github.com/virus-evolution/gofasta@latest && \
    cp /root/go/bin/gofasta /usr/local/bin/

# Instala o zip
RUN apt-get update && apt-get install -y zip

# instalar o gunicorn
RUN pip install gunicorn


# ---------- Instala dependências Python ----------
COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

# ---------- Copia o código ----------
COPY . .

# ---------- Expõe a porta da aplicação ----------
EXPOSE 5000

# ---------- Comando para iniciar a aplicação ----------
CMD ["python3", "-u", "app.py"]
