#!/usr/bin/env python
#!/usr/bin/env python
import numpy as np
import sys,pylab
sys.path.append('py/')
from DataStream import *
from BMXDemo import BMXDemoSingleDish
import cPickle as cP
telescope=BMXDemoSingleDish()
stream=DataStream(telescope)
## cut number of samples to 10 seconds
stream.tlist=stream.tlist[:10]
stream.fillStream()
cP.dump(stream,open("testStream.pickle",'w'))

