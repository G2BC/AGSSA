from flask import Flask
from celery import Celery
import os
from dotenv import load_dotenv

load_dotenv()

app = Flask(__file__)
celery = Celery()


aplicacao_path = os.getenv('APLICACAO_PATH')
AGSSA_treinamento_path = 'AGSSA.treinamento.AGSSA_treinamento'
AGSSA_analise_path = os.path.join(aplicacao_path, 'AGSSA/AGSSA_classificacao.py')


UPLOADS_PATH = os.getenv("UPLOADS_PATH", os.path.join(aplicacao_path, "uploads"))
RESULTS_PATH = os.getenv("RESULTS_PATH", os.path.join(aplicacao_path, "results"))