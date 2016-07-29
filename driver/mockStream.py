#!/usr/bin/env python
#!/usr/bin/env python
import numpy as np
import sys,pylab
sys.path.append('py/')
from DataStream import *
from BMXDemo import BMXDemoSingleDish

telescope=BMXDemoSingleDish()
stream=DataStream(telescope)
stream.fillStream()
