import celery
import traceback
from .CrimeReader import *
from .ObserveSky import *

app = celery.Celery('celery_tasks', broker="redis://astro0010:6379/0", backend="redis://astro0010:6379/1")
app.conf.accept_content = ['pickle']
app.conf.result_serializer='pickle'
app.conf.task_serializer='pickle'

logger = celery.utils.log.get_task_logger(__name__)

@app.task
def get_stream(telescope, skyc_list, nu, i, whichfield, psources):
    reader=CrimeReader()
    if whichfield is not None:
        # Read input map
        field=reader.named_slice(whichfield,i)
        # Generate time stream
        perfreqstreams=getIntegratedSignal(telescope, skyc_list, field, nu, Npix=201, Nfwhm=3)
    else:
        # Return Zero
        perfreqstreams=[np.zeros(len(skyc_list)) for i in range(len(telescope.beams))]

    # Add point sources if requested
    if (psources is not None):
        perfreqs=getPointSourceSignal(telescope, skyc_list, psources, nu)
        for i,s in enumerate(perfreqs):
            perfreqstreams[i]+=s

    return perfreqstreams
