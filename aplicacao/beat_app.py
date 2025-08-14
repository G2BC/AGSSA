from celery.schedules import schedule
from extensions import celery



celery.conf.beat_schedule = {
    'remove-old-directories-every-10-seconds': {
        'task': 'tasks.remove_many_directorys',
        'schedule': schedule(10.0),  
        'args': (),  
    },
}
