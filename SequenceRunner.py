import numpy as np
import CondSornB2SequenceDistRunner as CSR
import gc
import pdb

#CPR.clear(erase=True,all=True)
#CPR.reinit()
gc.collect()
n=5
ei=0.2
ie=0.2
delays=[0.000,0.001,0.002,0.003]
dists=[1,3,5,6]
#dist_num = 4
#dist_delay = 0.0 #in seconds
for run_num in range(n):
    for dist_num in dists:
    	for dist_delay in delays:
    		run_tag = "test_ei"+str(ei)+"_ie"+str(ie)+"D"+str(dist_num)+"del"+str(dist_delay)+"run"+str(run_num)
    		print(run_tag+":")
    		#pdb.set_trace()
    		CSR.run_it(run_tag,ei,ie,dist_num,dist_delay,play=False,save_it=True)
    		#CPR.clear(erase=True,all=True)
    		#CPR.reinit()
    		gc.collect()
gc.collect()
