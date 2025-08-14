from flask import Flask
from celery import Celery
import os
from dotenv import load_dotenv

load_dotenv()

app = Flask(__file__)
celery = Celery()



aplicacao_path = os.getenv('APLICACAO_PATH')
agua_treinamento_path = os.path.join(aplicacao_path, 'AGUA/AGUA_treinamento.py')
agua_analise_path = os.path.join(aplicacao_path, 'AGUA/AGUA_classificacao.py')


UPLOADS_PATH = os.getenv("UPLOADS_PATH", os.path.join(aplicacao_path, "uploads"))
RESULTS_PATH = os.getenv("RESULTS_PATH", os.path.join(aplicacao_path, "results"))