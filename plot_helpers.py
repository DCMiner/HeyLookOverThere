from brian2 import *
import peakutils as pu
from collections import Counter
from copy import deepcopy
from scipy import sparse
import analysis_helpers as AH
import matplotlib.pyplot as pt
import elephant.conversion as conv
import neo as ne
import quantities as pq
import pdb

def thresh_peak_times(peak_times,upper=1,lower=0):
	peaks_shape = shape(peak_times)
	threshed_peaks = []
	for i in range(peaks_shape[1]):
		temp_peaks_u = peak_times[:,i][peak_times[:,i]<upper]
		temp_peaks_ul = temp_peaks_u[temp_peaks_u >lower]
		threshed_peaks.append(temp_peaks_ul)
	return threshed_peaks

def peak_violin(pks,show=True): #assume 5 part linear sequence
	pts = pt.violinplot(pks,vert=False,showmeans=True,showextrema=False,points=50)
	pts['bodies'][0].set_facecolor('blue')
	pts['bodies'][0].set_edgecolor('blue')
	pts['bodies'][0].set_label('A')
	pts['bodies'][1].set_facecolor('orange')
	pts['bodies'][1].set_edgecolor('orange')
	pts['bodies'][1].set_label('B')
	pts['bodies'][2].set_facecolor('green')
	pts['bodies'][2].set_edgecolor('green')
	pts['bodies'][2].set_label('C')
	pts['bodies'][3].set_facecolor('red')
	pts['bodies'][3].set_edgecolor('red')
	pts['bodies'][3].set_label('D')
	pts['bodies'][4].set_facecolor('magenta')
	pts['bodies'][4].set_edgecolor('magenta')
	pts['bodies'][4].set_label('E')
	pt.yticks([1,2,3,4,5],['A','B','C','D','E'])
	if show:
		pt.xlim(0,0.008)
		pt.show()

def ctrl_peak_hist(peak_times,upper=1,legend=True,dens=True,show=True):
    pt.figure()
    pt.hist(peak_times[:,0][peak_times[:,0]<upper],bins=50,density=dens,label='A',histtype='step')
    pt.hist(peak_times[:,1][peak_times[:,1]<upper],bins=50,density=dens,label='B',histtype='step')
    pt.hist(peak_times[:,2][peak_times[:,2]<upper],bins=50,density=dens,label='C',histtype='step')
    pt.hist(peak_times[:,3][peak_times[:,3]<upper],bins=50,density=dens,label='D',histtype='step')
    pt.hist(peak_times[:,4][peak_times[:,4]<upper],bins=50,density=dens,label='E',histtype='step')
    if legend:
    	pt.legend()
    if show:
    	pt.xlabel('time (s)')
    	pt.ylabel('count (arbitrary units)')
    	pt.xlim(0.,0.007)
    	pt.tight_layout()
    	pt.show()


def peak_hist(dist,delay,peak_times_array,upper=1,legend=True,dens=True,show=True):
    pt.hist(peak_times_array[dist,delay,:,0][peak_times_array[dist,delay,:,0]<upper],bins=50,density=dens,label='A',histtype='step')
    pt.hist(peak_times_array[dist,delay,:,1][peak_times_array[dist,delay,:,1]<upper],bins=50,density=dens,label='B',histtype='step')
    pt.hist(peak_times_array[dist,delay,:,2][peak_times_array[dist,delay,:,2]<upper],bins=50,density=dens,label='C',histtype='step')
    pt.hist(peak_times_array[dist,delay,:,3][peak_times_array[dist,delay,:,3]<upper],bins=50,density=dens,label='D',histtype='step')
    pt.hist(peak_times_array[dist,delay,:,4][peak_times_array[dist,delay,:,4]<upper],bins=50,density=dens,label='E',histtype='step')
    pt.xlim(0.,0.007)
    if legend:
    	pt.legend()
    if show:
    	pt.show()

def peak_hists(dists,delays,peak_times_array,peak_vals_array,upper=1,legend=True,dens=True,show=True):
	dists_num = np.size(dists)
	dels_num =np.size(delays)
	success_rates = AH.get_success_rates(peak_times_array,peak_vals_array,thresh=10)
	counter = 0
	pt.figure() 
	for dist,delay in np.ndindex((dists_num,dels_num)):
		counter += 1
		pt.subplot(dists_num,dels_num,counter)
		peak_hist(dist,delay,peak_times_array,upper,False,dens,False)
		#title_string = "dist="+str(dist)+", del="+str(delay)
		##current trial structure specific:
		dist_str = ''
		if dist == 0:
			dist_str = 'A'
		elif dist == 1:
			dist_str = 'C'
		elif dist == 2:
			dist_str = 'E'
		elif dist == 3:
			dist_str = 'ext'
		passed = success_rates[dist,delay]
		title_string = "dist="+dist_str+", del="+str(delay)+" ms \n passed="+str(passed)
		###current trial structure specific:
		#pt.text(0.001,100,title_string)
		pt.text(0.0045,1000,title_string)
		pt.xlim(0.,0.007)
		if legend and counter == 1:
			pt.legend()
	if show:
		pt.tight_layout()
		pt.show()

def peak_hists_singledist(dist,delays,peak_times_array,peak_vals_array,upper=1,legend=True,dens=True,show=True):
	dists_num = 1
	dels_num =np.size(delays)
	success_rates = AH.get_success_rates(peak_times_array,peak_vals_array,thresh=10)
	counter = 0
	pt.figure() 
	for delay in range(dels_num):
		counter += 1
		pt.subplot(dels_num,dists_num,counter)
		peak_hist(dist,delay,peak_times_array,upper,False,dens,False)
		#title_string = "dist="+str(dist)+", del="+str(delay)
		##current trial structure specific:
		dist_str = ''
		if dist == 0:
			dist_str = 'A'
		elif dist == 1:
			dist_str = 'C'
		elif dist == 2:
			dist_str = 'E'
		elif dist == 3:
			dist_str = 'ext'
		passed = success_rates[dist,delay]
		title_string = "dist="+dist_str+", del="+str(delay)+" ms \n passed="+str(passed)
		###current trial structure specific:
		#pt.text(0.001,100,title_string)
		#pt.text(0.0056,150,title_string)
		pt.text(0.0045,1000,title_string)
		pt.xlim(0.,0.007)
		if legend and counter == 1:
			pt.legend()
	if show:
		pt.tight_layout()
		pt.show()

def plot_dis_dev(dists,delays,dis_mns,dis_vrs,dev_mns,dev_vrs):
	pt.figure()
	pt.subplot(221)
	pt.imshow(dis_mns)
	pt.title('disruption index mean')
	pt.xticks(np.arange(np.size(delays)),delays)
	pt.xlabel('delays (s)')
	pt.yticks(np.arange(np.size(dists)),dists)
	pt.ylabel('distractor position')
	pt.colorbar()
	pt.subplot(222)
	pt.imshow(dis_vrs)
	pt.title('disruption index variance')
	pt.xticks(np.arange(np.size(delays)),delays)
	pt.xlabel('delays (s)')
	pt.yticks(np.arange(np.size(dists)),dists)
	pt.ylabel('distractor position')
	pt.colorbar()
	pt.subplot(223)
	pt.imshow(dev_mns)
	pt.title('deviation index mean')
	pt.xticks(np.arange(np.size(delays)),delays)
	pt.xlabel('delays (s)')
	pt.yticks(np.arange(np.size(dists)),dists)
	pt.ylabel('distractor position')
	pt.colorbar()
	pt.subplot(224)
	pt.imshow(dev_vrs)
	pt.title('deviation index variance')
	pt.xticks(np.arange(np.size(delays)),delays)
	pt.xlabel('delays (s)')
	pt.yticks(np.arange(np.size(dists)),dists)
	pt.ylabel('distractor position')
	pt.colorbar()
	pt.tight_layout()
	pt.show()

def plot_dis_dev2hack(dists,delays,dis_mns,dev_mns):
	pt.figure()
	#pt.tight_layout()
	pt.subplot(121)
	pt.imshow(dis_mns)
	pt.title('disruption index mean')
	pt.xticks(np.arange(np.size(delays)),delays)
	#pt.xlabel('delay (ms)')
	#pt.yticks(np.arange(np.size(dists)),dists)
	pt.yticks(np.arange(np.size(dists)),['A','C','E','ext'])
	pt.ylabel('distractor position')
	pt.colorbar(fraction=0.046, pad=0.04)
	pt.subplot(122)
	pt.imshow(dev_mns)
	pt.title('deviation index mean')
	pt.xticks(np.arange(np.size(delays)),delays)
	pt.xlabel('distractor delay (s)')
	pt.yticks(np.arange(np.size(dists)),['A','C','E','ext'])
	#pt.ylabel('distractor position')
	pt.colorbar(fraction=0.046, pad=0.04)
	pt.tight_layout()
	pt.show()

def dist_pcadist(dists,delays,ptd_array,d_t=0.0005*second,max_t=0.05*second,legend=True,show=True):
	dists_num = np.size(dists)
	dels_num =np.size(delays)
	counter = 0
	all_steps = np.shape(ptd_array)[3]
	max_steps = int(max_t / d_t)
	lim_steps = all_steps
	if (all_steps > max_steps):
		lim_steps = max_steps
	lim_time = lim_steps * d_t
	time_steps = arange(0,lim_time,d_t)
	fig = pt.figure()

	ax = fig.add_subplot(111)
	ax.set_xlabel('time [seconds]')
	ax.set_ylabel('distance [arbitrary units]')
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

	for dist,delay in np.ndindex((dists_num,dels_num)):
		counter += 1
		fig.add_subplot(dists_num,dels_num,counter)
		pt.plot(time_steps,ptd_array[dist,delay].mean(0)[:lim_steps])
		#pt.plot(time_steps,ptd_array[delay,dist].mean(0)[:lim_steps])
		#title_string = "dist="+str(dist)+", del="+str(delay)
		##current trial structure specific:
		dist_str = ''
		if dist == 0:
			dist_str = 'A'
		elif dist == 1:
			dist_str = 'C'
		elif dist == 2:
			dist_str = 'E'
		elif dist == 3:
			dist_str = 'ext'
		title_string = "dist="+dist_str+", del="+str(delay)+" ms"
		###current trial structure specific:
		pt.text(0.002,100,title_string)
		#pt.text(0.002,1.0,title_string)
		#pt.ylim(0.,1.5)
		pt.ylim(0.,120)
		if legend and counter == 1:
			pt.legend()
	if show:
		pt.show()

def raster_scroll(binned_spike_counts,start=0,width=300,end=700000):
	import matplotlib.animation as animation
	fig = pt.figure()
	ims = []
	frames = end - start - width
	for i in range(frames):
		frame = pt.imshow(binned_spike_counts[:,i:i+width])
		ims.append(frame) #change time
	ani = animation.ArtistAnimation(fig, ims, interval=5, blit=True,repeat_delay=1000)
	pt.show()

def compare_peaks_dist_hist(peak_times_array,ctrl_times_array,test_index,del_index,legend=True,show=True):
	exp_shape = np.shape(peak_times_array)
	dists = exp_shape[0]
	dels = exp_shape[1]
	trials = exp_shape[2]
	stims = exp_shape[3]
	ctrl_shape = np.shape(ctrl_times_array)
	ctrl_trials = ctrl_shape[0]
	if (ctrl_shape[1] != stims) or (test_index > stims) or (del_index > dels):
		print('parameter mismatch')
		return
	exps = peak_times_array[:,del_index,:,:]
	bins = 50
	for i in range(dists):
		label = 'dist = ' + str(i)
		hist(exps[i,:,test_index][exps[i,:,test_index]<999],histtype='step',label=label,bins=bins)
	label = 'ctrl ' + str(test_index)
	hist(ctrl_times_array[:,test_index][ctrl_times_array[:,test_index]<999],histtype='step',label=label,bins=bins)
	if legend:
		pt.legend()
	if show:
		pt.tight_layout()
		pt.show()

def raster_plot(spike_times_list,start=0*second,end=350*second,dt=0.1*ms,units=True,show=True,shuffle=False):
	total_time = end - start
	steps = int(total_time / dt)
	neurons = len(spike_times_list.item())
	neuron_list = []
	spike_list = []
	indexlist = np.arange(0,neurons)
	if shuffle:
		np.random.shuffle(indexlist)
	indexcount = 0
	for i in indexlist:
		for j in range(len(spike_times_list.item()[i])):
			spike = spike_times_list.item()[i][j]
			if (spike >= start) and (spike <= end):
				#neuron_list.append(i)
				neuron_list.append(indexcount)
				spike_list.append(spike)
		indexcount += 1
	pt.plot(spike_list,neuron_list,'.')
	pt.xlim(start/second,end/second)
	ylabel('neuron number')
	xlabel('time (s)')
	if show:
		pt.tight_layout()
		pt.show()

def raster_plot_colorhack(spike_times_list,start=0*second,end=350*second,dt=0.1*ms,units=True,show=True):
	total_time = end - start
	steps = int(total_time / dt)
	neurons = len(spike_times_list.item())
	neuron_list = []
	spike_list = []
	for i in range(neurons):
		for j in range(len(spike_times_list.item()[i])):
			spike = spike_times_list.item()[i][j]
			if (spike >= start) and (spike <= end):
				neuron_list.append(i)
				spike_list.append(spike)
	color_list = []
	for n in neuron_list:
		if (n >= 0) and (n < 20):
			color_list.append('blue')
		elif (n >= 20) and (n < 40):
			color_list.append('orange')
		elif (n >= 40) and (n < 60):
			color_list.append('green')
		elif (n >= 60) and (n < 80):
			color_list.append('red')
		elif (n >= 80) and (n < 100):
			color_list.append('magenta')
		else:
			color_list.append('black')	
	fig, ax1 = pt.subplots()
	ax1.scatter(spike_list,neuron_list,color=color_list,s=2.)
	ax1.set_xlim(start/second,end/second)
	ax1.set_xlabel('time (s)')
	ax1.set_ylabel('neuron number')
	ax2 = ax1.twinx()
	ax2.set_ylabel('neuron group')
	#ax2.tick_params('y', colors='b')
	ax2.set_yticks([0.05,0.15,0.25,0.35,0.45,0.75])
	ax2.set_yticklabels(['A','B','C','D','E','ext'])
	fig.tight_layout()
	if show:
		pt.show()

def compareshufflerasters(spike_times_list,start=199.998*second,end=200.008*second,dt=0.1*ms,units=True):
	total_time = end - start
	steps = int(total_time / dt)
	neurons = len(spike_times_list.item())
	neuron_list = []
	spike_list = []
	indexlist = np.arange(0,neurons)
	if shuffle:
		np.random.shuffle(indexlist)
	indexcount = 0
	for i in indexlist:
		for j in range(len(spike_times_list.item()[i])):
			spike = spike_times_list.item()[i][j]
			if (spike >= start) and (spike <= end):
				#neuron_list.append(i)
				neuron_list.append(indexcount)
				spike_list.append(spike)
		indexcount += 1
	pt.subplot(1,2,1)
	pt.plot(spike_list,neuron_list,'.')
	pt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
	pt.xlim(start/second,end/second)
	pt.xlabel('time (s)')
	pt.ylabel('neuron number (unsorted)')
	neuron_list = []
	spike_list = []
	for i in range(neurons):
		for j in range(len(spike_times_list.item()[i])):
			spike = spike_times_list.item()[i][j]
			if (spike >= start) and (spike <= end):
				neuron_list.append(i)
				spike_list.append(spike)
	color_list = []
	for n in neuron_list:
		if (n >= 0) and (n < 20):
			color_list.append('blue')
		elif (n >= 20) and (n < 40):
			color_list.append('orange')
		elif (n >= 40) and (n < 60):
			color_list.append('green')
		elif (n >= 60) and (n < 80):
			color_list.append('red')
		elif (n >= 80) and (n < 100):
			color_list.append('magenta')
		else:
			color_list.append('black')
	pt.subplot(1,2,2)
	pt.scatter(spike_list,neuron_list,color=color_list,s=2.)
	pt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
	pt.xlim(start/second,end/second)
	pt.xlabel('time (s)')
	pt.ylabel('neuron number (sorted)')
	pt.tight_layout()
	pt.show()



def peak_hists_singledist_withraster(dist,delays,peak_times_array,peak_vals_array,spike_times,trial_time,trial_length=7*ms,upper=1,legend=True,dens=True,show=True):
	#dists_num = 1
	dels_num =np.size(delays)
	del_nums = (np.array(delays)*1000).astype('int') #hacky  
	success_rates = AH.get_success_rates(peak_times_array,peak_vals_array,thresh=10)
	counter = 0
	pt.figure() 
	#for delay in range(dels_num):
	for delay in del_nums:
		counter += 1
		pt.subplot(dels_num,2,counter)
		peak_hist(dist,delay,peak_times_array,upper,False,dens,False)
		#title_string = "dist="+str(dist)+", del="+str(delay)
		##current trial structure specific:
		dist_str = ''
		if dist == 0:
			dist_str = 'A'
		elif dist == 1:
			dist_str = 'C'
		elif dist == 2:
			dist_str = 'E'
		elif dist == 3:
			dist_str = 'ext'
		passed = success_rates[dist,delay]
		delstr = str(delay)
		title_string = "dist="+dist_str+", del="+delstr+" ms \n passed="+str(passed)
		###current trial structure specific:
		#pt.text(0.001,100,title_string)
		#pt.text(0.0056,150,title_string)
		pt.xlim(0.,0.007)
		if dels_num == 2:
			pt.text(0.001,3000,title_string)
			pt.ylim(0,4000)
			if counter == 3:
				pt.xlabel('time (s)')
				pt.ylabel('count (arbitrary units)')
				if legend:
					pt.legend()
		else:
			pt.text(0.0045,1000,title_string)
			if counter == 7:
				pt.xlabel('time (s)')
				pt.ylabel('count (arbitrary units)')
				if legend:
					pt.legend()
		counter += 1
		pt.subplot(dels_num,2,counter)
		end_time = trial_time + trial_length
		dt = 0.1 * ms #hacky
		run = 0 #hacky
		steps = int(trial_time / dt)
		neurons = len(spike_times[dist][delay][run].item())
		neuron_list = []
		spike_list = []
		for i in range(neurons):
			for j in range(len(spike_times[dist][delay][run].item()[i])):
				spike = spike_times[dist][delay][run].item()[i][j]
				if (spike >= trial_time) and (spike <= end_time):
					neuron_list.append(i)
					spike_list.append(spike)
		color_list = []
		for n in neuron_list:
			if (n >= 0) and (n < 20):
				color_list.append('blue')
			elif (n >= 20) and (n < 40):
				color_list.append('orange')
			elif (n >= 40) and (n < 60):
				color_list.append('green')
			elif (n >= 60) and (n < 80):
				color_list.append('red')
			elif (n >= 80) and (n < 100):
				color_list.append('magenta')
			else:
				color_list.append('black')	
		pt.scatter(spike_list,neuron_list,color=color_list,s=2.)
		pt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
		pt.xlim(trial_time/second,end_time/second)
		if dels_num == 2:
			if counter == 4:
				pt.xlabel('time (s)')
				pt.ylabel('neuron number')
		else:
			if counter == 8:
				pt.xlabel('time (s)')
				pt.ylabel('neuron number')

	if show:
		pt.tight_layout()
		pt.show()

def peak_vio_singledist_withraster(dist,delays,peak_times_array,peak_vals_array,spike_times,trial_time,trial_length=7*ms,upper=1,legend=True,dens=True,show=True):
	#dists_num = 1
	dels_num =np.size(delays)
	del_nums = (np.array(delays)*1000).astype('int') #hacky  
	success_rates = AH.get_success_rates(peak_times_array,peak_vals_array,thresh=10)
	counter = 0
	pt.figure() 
	#for delay in range(dels_num):
	for delay in del_nums:
		counter += 1
		pt.subplot(dels_num,2,counter)
		#peak_hist(dist,delay,peak_times_array,upper,False,dens,False)
		pks = thresh_peak_times(peak_times_array[dist,delay])
		peak_violin(pks,show=False)
		#title_string = "dist="+str(dist)+", del="+str(delay)
		##current trial structure specific:
		dist_str = ''
		if dist == 0:
			dist_str = 'A'
		elif dist == 1:
			dist_str = 'C'
		elif dist == 2:
			dist_str = 'E'
		elif dist == 3:
			dist_str = 'ext'
		passed = success_rates[dist,delay]
		delstr = str(delay)
		passed_str = "passed=%s"%float('%.3g'%(passed*100))+"%"
		#title_string = "dist="+dist_str+", del="+delstr+" ms \n"+passed_str
		title_string = passed_str
		###current trial structure specific:
		#pt.text(0.001,100,title_string)
		#pt.text(0.0056,150,title_string)
		pt.xlim(0.,0.007)
		if dels_num == 2:
			pt.text(0.003,0.75,title_string)
			#pt.ylim(0,4000)
			if counter == 1:
				pt.xticks([])
				pt.ylabel('delay='+delstr+' ms \n distractor='+dist_str)
			if counter == 3:
				pt.xlabel('time (s)')
				pt.ylabel('delay='+delstr+' ms \n neuron group')
				#if legend:
					#pt.legend()
		else:
			pt.text(0.003,0.75,title_string)
			if counter == 7:
				pt.xlabel('time (s)')
				pt.ylabel('neuron group')
				#if legend:
					#pt.legend()
		counter += 1
		pt.subplot(dels_num,2,counter)
		end_time = trial_time + trial_length
		dt = 0.1 * ms #hacky
		run = 0 #hacky
		steps = int(trial_time / dt)
		neurons = len(spike_times[dist][delay][run].item())
		neuron_list = []
		spike_list = []
		for i in range(neurons):
			for j in range(len(spike_times[dist][delay][run].item()[i])):
				spike = spike_times[dist][delay][run].item()[i][j]
				if (spike >= trial_time) and (spike <= end_time):
					neuron_list.append(i)
					spike_list.append(spike)
		color_list = []
		for n in neuron_list:
			if (n >= 0) and (n < 20):
				color_list.append('blue')
			elif (n >= 20) and (n < 40):
				color_list.append('orange')
			elif (n >= 40) and (n < 60):
				color_list.append('green')
			elif (n >= 60) and (n < 80):
				color_list.append('red')
			elif (n >= 80) and (n < 100):
				color_list.append('magenta')
			else:
				color_list.append('black')	
		pt.scatter(spike_list,neuron_list,color=color_list,s=2.)
		pt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
		pt.gca().get_yaxis().tick_right()
		pt.xlim(trial_time/second,end_time/second)
		if dels_num == 2:
			if counter == 2:
				pt.xticks([])
			if counter == 4:
				pt.xlabel('time (s)')
				pt.ylabel('neuron number')
		else:
			if counter == 8:
				pt.xlabel('time (s)')
				pt.ylabel('neuron number')

	pt.tight_layout()
	if show:
		pt.show()

def ctrl_hist_vio(peak_times):
	pks = thresh_peak_times(peak_times)
	pt.figure()
	pt.subplot(121)
	peak_violin(pks,show=False)
	pt.xlim(0.,0.007)
	pt.xlabel('time (s)')
	pt.ylabel('neuron group')
	pt.subplot(122)
	#ctrl_peak_hist(peak_times,legend=True,dens=True,show=True)
	pt.hist(pks[0],bins=10,density=True,label='A',histtype='step')
	pt.hist(pks[1],bins=10,density=True,label='B',histtype='step')
	pt.hist(pks[2],bins=10,density=True,label='C',histtype='step')
	pt.hist(pks[3],bins=10,density=True,label='D',histtype='step')
	pt.hist(pks[4],bins=10,density=True,label='E',histtype='step')
	pt.legend()
	#pt.xlabel('time (s)')
	pt.ylabel('probability density')
	pt.xlim(0.,0.007)
	pt.tight_layout()
	pt.show()

def all_dist_vio_plots(peak_times,peak_vals,spike_times,trial_time=200*second):
	peak_vio_singledist_withraster(0,[0.000,0.001],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(0,[0.002,0.003],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)
	peak_vio_singledist_withraster(1,[0.000,0.001],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(1,[0.002,0.003],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(2,[0.000,0.001],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(2,[0.002,0.003],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(3,[0.000,0.001],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)	
	peak_vio_singledist_withraster(3,[0.002,0.003],peak_times,peak_vals,spike_times,trial_time,trial_length=7*ms,show=False)
	pt.show()

def rep_weights_hist(weight_stack,clust_size=20,active_clust_n=5,normed=True,show=True):
	aw,rw,of,nf,ob,nb,te,fe = AH.rep_weights_from_stack(weight_stack,clust_size,active_clust_n)
	pt.figure()
	pt.subplot(231)
	of = np.array(of) * 1e9
	pt.hist(of,density=normed)
	tmp_mean = np.mean(of)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(of)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'one forward\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.subplot(232)
	nf = np.array(nf) * 1e9
	pt.hist(nf,density=normed)
	tmp_mean = np.mean(nf)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(nf)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'n forward\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.subplot(233)
	te = np.array(te) * 1e9
	pt.hist(te,density=normed)
	tmp_mean = np.mean(te)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(te)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'to external\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.subplot(234)
	ob = np.array(ob) * 1e9
	pt.hist(ob,density=normed)
	tmp_mean = np.mean(ob)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(ob)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'one backward\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	if normed:
		pt.ylabel('probability density')
	else:
		pt.ylabel('count (arbitrary units)')
	pt.subplot(235)
	nb = np.array(nb) * 1e9
	pt.hist(nb,density=normed)
	tmp_mean = np.mean(nb)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(nb)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'n backward\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.subplot(236)
	fe = np.array(fe) * 1e9
	pt.hist(fe,density=normed)
	tmp_mean = np.mean(fe)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(fe)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'from external\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.xlabel('weight (nS)')
	pt.tight_layout()

	pt.figure()
	pt.subplot(121)
	aw = np.array(aw) * 1e9
	pt.hist(aw,density=normed)
	tmp_mean = np.mean(aw)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(aw)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'all weights\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	if normed:
		pt.ylabel('probability density')
	else:
		pt.ylabel('count (arbitrary units)')
	pt.subplot(122)
	rw = np.array(rw) * 1e9
	pt.hist(rw,density=normed)
	tmp_mean = np.mean(rw)
	mean_str = 'mean = %.3f nS' %tmp_mean
	tmp_std =  np.std(rw)
	std_str = 'std = %.3f nS' %tmp_std
	title_string = 'recurrent weights\n' + mean_str + '\n' + std_str
	pt.title(title_string)
	pt.xlabel('weight (nS)')

	if show:
		pt.tight_layout()
		pt.show()
