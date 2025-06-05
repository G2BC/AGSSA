# Trabalho de Conclusão de Curso - Joílson

## Table of Contents

- [Sobre](#sobre)
- [Começando](#começando)
- [Usando](#usando)

## [Sobre](#sobre)

Repositório para o desenvolvimento do Trabalho de Conclusão de Curso do curso de Análise e Desenvolvimento de Sistemas da Universidade do Estado da Bahia.

## [Começando](começando)

Instruções para instalação e execução do projeto.

### Pré-requisitos

O que você precisa para instalar o software e como instalá-lo.

- Python 3.8
- RabbitMQ (sudo apt install rabbitmq-server)
- Alinhamento das Sequências
  - Golang (<https://go.dev/doc/install>)
  - Minimap2 (<https://github.com/lh3/minimap2?tab=readme-ov-file#installation>)
  - Gofasta (<https://github.com/virus-evolution/gofasta?tab=readme-ov-file#installation>)

### Executando

Criando o ambiente virtual:

```bash
python3 -m venv venv
```

Ativando o ambiente virtual:

```bash
source venv/bin/activate
```

Instalando as dependências:


```bash
pip install -r aplicacao/requirements.txt
```

Criar arquivo .env com as infos que estão no arquivo .env.example

Executando o celery:

```bash
celery -A app.celery worker --loglevel=info
```

Executando a aplicação para desenvolvimento:

```bash
flask --app app run --debug
```

## [Usando](usando)

Ao executar o projeto o mesmo será executado em <http://localhost:5000>
