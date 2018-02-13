import celery
import traceback
from MapReader import *
from ObserveSky import *
from DataStream import get_stream

app = celery.Celery('celery_tasks', broker="redis://astro0010:6379/0", backend="redis://astro0010:6379/1")
app.conf.accept_content = ['pickle']
app.conf.result_serializer='pickle'
app.conf.task_serializer='pickle'

logger = celery.utils.log.get_task_logger(__name__)

@app.task
def get_stream_celery(reader, telescope, skyc_list, i, fld, psources, Npix, Nfwhm):
    return get_stream(reader, telescope, skyc_list, i, fld, psources, Npix, Nfwhm)
