from flask import Flask
from celery import Celery


app = Flask(__file__)
celery = Celery()

