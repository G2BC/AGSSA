from celery.schedules import crontab
from extensions import celery
import tasks  

celery.conf.timezone = 'America/Sao_Paulo'
celery.conf.enable_utc = False  

celery.conf.beat_schedule = {
    'remove-old-directories': {
        'task': 'tasks.remove_old_directorys',
        'schedule': crontab(hour=23, minute=59),  
        'args': (),
    },
}
