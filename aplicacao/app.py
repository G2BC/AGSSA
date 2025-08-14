import os
from dotenv import load_dotenv
from extensions import app, celery
import tasks  
import routes

load_dotenv()

# Config Flask
app.secret_key = os.getenv('SECRET_KEY')
app.config['CELERY_BROKER_URL'] = os.getenv('CELERY_BROKER_URL')
app.config['CELERY_RESULT_BACKEND'] = os.getenv('CELERY_RESULT_BACKEND')

# Config Celery
celery.conf.update(app.config)
celery.conf.broker_url = app.config['CELERY_BROKER_URL']
celery.conf.result_backend = app.config['CELERY_RESULT_BACKEND']


if __name__ == '__main__':
    app.run(debug=True)
