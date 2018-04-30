from datamanager import datamanager
import farmit
import numpy as np

dm = datamanager()
dm.gettags(reduced=True)
year,month,day,hr,min = dm.parsetags(dm.tags)
ind = np.where((year==18) & (month==4) & (day>=1) & (day<=2))[0]
tags = dm.tags[ind]

fields = ['colore','egfree','gfree','gsync','psources','colore+egfree+gfree+gsync+psources']
decdither = 0.5

sn = ['00000', '00001']
nu_f = [800.0, np.inf]

for k in range(len(sn)):
    for tag in tags:
        f = farmit.farmit('driver/mockStream.py', 
                          args={'tag':[tag], 'field':fields, 'sn':[sn[k]], 'nu_f':[nu_f[k]],
                                'decdither':[decdither]},
                          OMP_NUM_THREADS=1)
        f.writejobfiles()


