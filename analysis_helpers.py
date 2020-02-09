from brian2 import *
import peakutils as pu
from collections import Counter
from copy import deepcopy
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import matplotlib.pyplot as pt
import elephant.conversion as conv
import neo as ne
import quantities as pq
import gc
import pdb

def get_rate(monitor, width=2*ms, window='gaussian'):
#extracts rate from monitor
	return monitor.smooth_rate(width=width,window=window)

def get_peaks(rate, thresh=0.025,subsamp=10,min_dist=1000): #extracts peaks from rate, converts time
	min_dist /= subsamp
	#rate = rate / Hz
	peaks = pu.indexes(rate[::subsamp], thres=thresh, min_dist=min_dist)
	peaks = peaks / (1e4/subsamp)
	return peaks.tolist()

#def truncate_peaks(peak_list,min_time=0, max_time=total_time):
def truncate_peaks(peak_list,min_time, max_time):
#gets truncated time of single peak list
	temp_peak_list = deepcopy(peak_list)
	while temp_peak_list[0] < min_time:
		temp_peak_list.pop(0)
	while temp_peak_list[len(temp_peak_list)-1] > max_time:
		temp_peak_list.pop(len(temp_peak_list)-1)
	return temp_peak_list

#def get_peak_list(mon_list,min_time=total_time-test_time, max_time=total_time):
def get_peak_list(mon_list,min_time,max_time):
#gets list of peak lists
	peak_list = []
	for i in range(len(mon_list)):
		peaks = truncate_peaks(get_peaks(get_rate(mon_list[i])),min_time,max_time)
		#peaks = get_peaks(get_rate(mon_list[i]))
		peak_list.append(peaks)
	return peak_list

def get_min_first_from_peak_list(peak_list):
#looks at bottom of list of peak lists, returns index of earliest
	temp_first = zeros(len(peak_list))
	for i in range(len(peak_list)):
		temp_first[i] = peak_list[i][0]
	return argmin(temp_first)

def check_empty(peak_list,offset=0):
#makes sure all lists in list of peak lists have entries
	empt = False
	for i in range(len(peak_list)):
		if len(peak_list[i]) < (1 + offset):
			empt = True
	return empt

def norm_mat(mat,axis=0):
    return mat / mat.sum(axis)

def make_trans_matrix(peak_list):
#run over list of peak lists, make transition matrix
	n_c = shape(peak_list)[0]
	trans_count_matrix = zeros(shape=(n_c,n_c))
	norm_count_matrix = zeros(shape=(n_c,n_c))
	temp_list = deepcopy(peak_list)
	#temp_list = peak_list
	current = get_min_first_from_peak_list(temp_list)
	while not check_empty(temp_list,1):
		temp_list[current].pop(0)
		nxt = get_min_first_from_peak_list(temp_list)
		norm_count_matrix[current] += ones(shape=(n_c))
		trans_count_matrix[current,nxt] += 1.0
		current = nxt
	#return trans_count_matrix / norm_count_matrix
	return norm_mat(trans_count_matrix / norm_count_matrix)

def monsmat(k,syns,monitor,num,sprs=True):
    temp_mat = sparse.lil_matrix((num,num))
    for i,j,w in zip(syns.i,syns.j,monitor[:,k]):
        temp_mat[i,j] = w / siemens
    if sprs: return temp_mat
    else: return temp_mat.todense()

def mean_ff(w_mat,rel=True,el=5,oc=0.5):
#w_mat is recurrent / square matrix
#returns mean feedforward weight
#if rel, relative to mean weight
#assumes 5 element sequence oppcupying first half of network
	n = shape(w_mat)[0]
	n_in = int (n*oc/el)
	w_mean = mean(w_mat[w_mat>0])
	w_ff_list = []
	for i in range(el-1):
		w_ff_temp = w_mat[n_in*i:n_in*(i+1),n_in*(i+1):n_in*(i+2)]
		w_ff_temp_list = w_ff_temp[w_ff_temp>0]
		for w in w_ff_temp_list.T:
			w_ff_list.append(w.item())
	w_ff_mean = mean(w_ff_list)
	if rel:
		return w_ff_mean / w_mean
	else:
		return w_ff_mean

#syns=syn_ee,monitor=weights_ee.w_ampa,num=n_e by default
def mean_ff_array(syns,monitor,num,rel=True):
	#print('a')
	n_t = shape(monitor)[1]
	#print('b')
	mean_ff_list = []
	#print('c')
	for i in range(n_t):
		#print('d1')
		w_temp = array(monsmat(i,syns,monitor,num,sprs=False))
		#print(w_temp)
		#print('d2')
		mean_ff_temp = mean_ff(w_temp,rel=rel)
		#print('d3')
		mean_ff_list.append(mean_ff_temp)
		#print('d4')
	#print('e')
	return array(mean_ff_list)

###assumes mon list ABCD for decisions ABC/ABD, min and max time are cue (A or partial A) triggered range
#def simple_decision_readout(mon_list,min_time,max_time,dbg=False):
def simple_decision_readout(rate_list,min_time,max_time,dbg=False):
	peak_list = []
	for i in range(4):
		peaks = truncate_peaks(get_peaks(rate_list[i],subsamp=10,min_dist=250),min_time,max_time)
		peak_list.append(peaks)
	peak_vals = []
	for i in range(len(peak_list)):
		temp_peak_vals = []
		for j in range(len(peak_list[i])):
			rate_index = int(peak_list[i][j]*10000) #assumes 0.1 ms dt
			temp_peak_vals.append(rate_list[i][rate_index])
		peak_vals.append(temp_peak_vals)
	#return peak_list,peak_vals #debug return
	#thresh values for "counting" ABCD
	thresh_list= [50.0,10.5,10.5,10.5]
	threshed_times = []
	threshed_indices = []
	for i in range(len(peak_list)):
		temp_threshed_times = []
		temp_threshed_indices = []	
		for j in range(len(peak_list[i])):
			#if (peak_vals[i][j] / hertz)>thresh_list[i]:
			if peak_vals[i][j]>thresh_list[i]:
				temp_threshed_times.append(peak_list[i][j])
				temp_threshed_indices.append(j)
		threshed_times.append(temp_threshed_times)
		threshed_indices.append(temp_threshed_indices)
	#return threshed_times,threshed_indices
	decisions = []
	for i in range(len(threshed_indices[0])):
		#decisions u-undecided, nb-no b, ncd-no decision, c-C, d-D
		decision = 'u' #P
		t_a = threshed_times[0][i]
		t_b_before_mask = array(threshed_times[1]) < t_a
		t_b_next_index = sum(t_b_before_mask*1)
		if t_b_next_index >= len(threshed_times[1]):
			t_b = max_time + 1.0
		else:
			t_b = threshed_times[1][t_b_next_index]
		if ((t_b - t_a) > 0.025):
			decision = 'nb'
		if (decision == 'u'):
			t_c_before_mask = array(threshed_times[2]) < t_b
			t_c_next_index = sum(t_c_before_mask*1)
			if t_c_next_index >= len(threshed_times[2]):
				t_c = max_time + 1.0
			else:
				t_c = threshed_times[2][t_c_next_index]
			t_d_before_mask = array(threshed_times[3]) < t_b
			t_d_next_index = sum(t_d_before_mask*1)
			if t_d_next_index >= len(threshed_times[3]):
				t_d = max_time + 1.0
			else:
				t_d = threshed_times[3][t_d_next_index]
			if ((t_c - t_b) > 0.025) and ((t_d - t_b) > 0.025):
				decision = 'ncd'
			else:
				if (t_d > max_time):
					decision = 'c'
				elif (t_c > max_time):
					decision = 'd'
				else:
					val_c = peak_vals[2][threshed_indices[2][t_c_next_index]]
					val_d = peak_vals[3][threshed_indices[3][t_d_next_index]]
					if (val_c > val_d):
						decision = 'c'
					else:
						decision = 'd'
		decisions.append(decision)
	if dbg:
		return peak_list,peak_vals,threshed_times,threshed_indices,decisions
	else:
		return decisions

def counts_from_dec_stack_entry(current_stack):
	count_stack = []
	c_stack = []
	d_stack = []
	nb_stack = []
	ncd_stack = []
	for i in range(len(current_stack)):
		count = Counter(current_stack[i])
		count_stack.append(count)
		c_stack.append(count['c'])
		d_stack.append(count['d'])
		nb_stack.append(count['nb'])
		ncd_stack.append(count['ncd'])
	c_stack = np.array(c_stack)
	d_stack = np.array(d_stack)
	nb_stack = np.array(nb_stack)
	ncd_stack = np.array(ncd_stack)
	return c_stack,d_stack,nb_stack,ncd_stack

def dec_arrays(dec_stack):
	coverd_mean = np.zeros((len(dec_stack),len(dec_stack[0])))
	coverd_var = np.zeros((len(dec_stack),len(dec_stack[0])))
	c_mean = np.zeros((len(dec_stack),len(dec_stack[0])))
	c_var = np.zeros((len(dec_stack),len(dec_stack[0])))
	d_mean = np.zeros((len(dec_stack),len(dec_stack[0])))
	d_var = np.zeros((len(dec_stack),len(dec_stack[0])))
	nb_mean = np.zeros((len(dec_stack),len(dec_stack[0])))
	nb_var = np.zeros((len(dec_stack),len(dec_stack[0])))
	ncd_mean = np.zeros((len(dec_stack),len(dec_stack[0])))
	ncd_var = np.zeros((len(dec_stack),len(dec_stack[0])))
	for i_ei in range(len(dec_stack)):
		for i_ie in range(len(dec_stack[0])):
			c_s,d_s,nb_s,ncd_s = counts_from_dec_stack_entry(dec_stack[i_ei,i_ie])
			coverd_mean[i_ei,i_ie] = (c_s / d_s).mean()
			coverd_var[i_ei,i_ie] = (c_s / d_s).std()	
			c_mean[i_ei,i_ie] = c_s.mean()
			c_var[i_ei,i_ie] = c_s.std()
			d_mean[i_ei,i_ie] = d_s.mean()
			d_var[i_ei,i_ie] = d_s.std()
			nb_mean[i_ei,i_ie] = nb_s.mean()
			nb_var[i_ei,i_ie] = nb_s.std()
			ncd_mean[i_ei,i_ie] = ncd_s.mean()
			ncd_var[i_ei,i_ie] = ncd_s.std()
	return coverd_mean,coverd_var,c_mean,c_var,d_mean,d_var,nb_mean,nb_var,ncd_mean,ncd_var

#new
def trial_rates(trial_rl,min_time=150.,max_time=250.,cue_freq=2.,window=0.01,active=5,dbg=False):
	dt = 10000
	runs = np.shape(trial_rl)[0]
	trials = int((max_time - min_time) * cue_freq)
	starts = np.arange(min_time,max_time,(1 / cue_freq))
	stops = np.arange(min_time + window,max_time + window,(1 / cue_freq))
	trial_traces = []
	for cue in range(trials):
		start = starts[cue]
		stop = stops[cue]
		for run in range(runs):
			traces = trial_rl[run,:,int(start*dt):int(stop*dt)]
			trial_traces.append(traces)
	return np.array(trial_traces)

#new
def sequential_peak_readout(trial_rl,min_time=150.,max_time=250.,cue_freq=2.,window=0.035,offset=-0.01,
	active=5):
	dt = 10000
	runs = np.shape(trial_rl)[0]
	trials = int((max_time - min_time) * cue_freq)
	starts = np.arange(min_time + offset,max_time + offset,(1 / cue_freq))
	stops = np.arange(min_time + offset + window,max_time + offset + window,(1 / cue_freq))
	peak_time_list = zeros((trials * runs,active))
	peak_val_list = zeros((trials * runs,active))
	t_count = 0
	for cue in range(trials):
		start = starts[cue]
		stop = stops[cue]
		for run in range(runs):
			traces = trial_rl[run,:,int(start*dt):int(stop*dt)]
			for cluster in range(active):
				trace = traces[cluster]
				peak_pos = get_peaks(trace,thresh=0,subsamp=1,min_dist=0)
				if np.size(peak_pos) < 1:
					peak_pos = 999 #failure flag
					peak_val = 0 #failute flag
				elif np.size(peak_pos) > 1:
					mult_peak_vals = trace[np.array(peak_pos * dt,dtype=int)]
					peak_index = np.argmax(mult_peak_vals)
					peak_pos = peak_pos[peak_index]
					peak_val = trace[int(peak_pos * dt)]
				else:
					peak_pos = peak_pos[0]
					peak_val = trace[int(peak_pos * dt)]
				#peak_time_list[t_count,cluster] = peak_pos
				peak_time_list[t_count,cluster] = peak_pos + offset 
				peak_val_list[t_count,cluster] = peak_val
			t_count += 1
	return peak_time_list,peak_val_list

#new
def sequence_list_veto(peak_time_list,peak_val_list,occ_check=True,thresh=0):
	vetod_time_list = []
	vetod_val_list = []
	for i in range(np.shape(peak_time_list)[0]):
		if occ_check and (sum(peak_time_list[i] > 1) == 0):
			vetod_time_list.append(peak_time_list[i])
			vetod_val_list.append(peak_val_list[i])
		elif not occ_check:
			vetod_time_list.append(peak_time_list[i])
			vetod_val_list.append(peak_val_list[i])
	vetod_time_list = np.array(vetod_time_list)
	vetod_val_list = np.array(vetod_val_list)
	threshed_time_list = []
	threshed_val_list = []
	for i in range(np.shape(vetod_time_list)[0]):
		if sum(vetod_val_list[i] > thresh) == np.shape(vetod_val_list)[1]:
			threshed_time_list.append(vetod_time_list[i])
			threshed_val_list.append(vetod_val_list[i])
	return np.array(threshed_time_list),np.array(threshed_val_list)

#new
def get_success_rate(peak_times,peak_vals,occ_check=True,thresh=0):
	trials = np.shape(peak_vals)[0]
	output = sequence_list_veto(peak_times,peak_vals,occ_check,thresh)
	successes = np.shape(output)[1]
	success_rates = successes / trials
	return success_rates

def get_success_rates(peak_time_array,peak_val_array,occ_check=True,thresh=0):
	array_shape = np.shape(peak_val_array)[0:2]
	trials = np.shape(peak_val_array)[2]
	clusters = np.shape(peak_val_array)[3]
	success_rates = np.zeros((array_shape))
	for i,j in np.ndindex(array_shape):
		success_rates[i,j] = get_success_rate(peak_time_array[i,j],peak_val_array[i,j],occ_check,thresh)
	return success_rates


#new
def get_mean_times(peak_time_list,peak_val_list,thresh=10,occ_check=True,return_vals=False):
	good_time_list,good_val_list = sequence_list_veto(peak_time_list,peak_val_list,occ_check,thresh)
	if return_vals:
		return good_time_list.mean(0),good_time_list.var(0),good_val_list.mean(0),good_val_list.var(0)
	else:
		return good_time_list.mean(0),good_time_list.var(0)

#new
def time_diffs_array(time_list):
	diffs = np.zeros((np.size(time_list)-1))
	for i in range(np.size(diffs)):
		diffs[i] = time_list[i+1] - time_list[i]
	return diffs

#new
def get_diffs(peak_time_list,peak_val_list,occ_check=True,thresh=10):
	good_time_list,good_val_list = sequence_list_veto(peak_time_list,peak_val_list,occ_check,thresh)
	diff_list = np.zeros((np.shape(good_time_list)[0],np.shape(good_time_list)[1]-1))
	for i in range(np.shape(good_time_list)[0]):
		diff_list[i] = time_diffs_array(good_time_list[i])
	return diff_list

#new
#invented measure
def disruption_index(ctrl_diffs_mean,ctrl_diffs_var,exp_diffs):
	mean_measure = (exp_diffs - ctrl_diffs_mean) / np.sqrt(ctrl_diffs_var)
	#disruption_index = sum(abs(mean_measure)) / np.size(ctrl_diffs_mean)
	disruption_index = sum(mean_measure) / np.size(ctrl_diffs_mean)
	return disruption_index

#new
#single condition experimental list
def mean_disruption_index(ctrl_diffs_mean,ctrl_diffs_var,exp_peak_times,exp_peak_vals,occ_check=True,
	thresh=10):
	disruption_index_mean = zeros((np.size(ctrl_diffs_mean)))
	disruption_index_var = zeros((np.size(ctrl_diffs_var)))
	diffs = get_diffs(exp_peak_times,exp_peak_vals,occ_check,thresh)
	disruption_index_list = zeros(np.shape(diffs)[0])
	for i in range(np.shape(diffs)[0]):
		disruption_index_list[i] = disruption_index(ctrl_diffs_mean,ctrl_diffs_var,diffs[i])
	return disruption_index_list.mean(),disruption_index_list.var()

#new
def scan_disruption(ctrl_peak_times,ctrl_peak_vals,exp_times_array,exp_vals_array,occ_check=True,
	thresh=10):
	scan_mean_matrix = np.zeros(np.shape(exp_times_array)[0:2])
	scan_var_matrix = np.zeros(np.shape(exp_times_array)[0:2])
	ctrl_diffs = get_diffs(ctrl_peak_times,ctrl_peak_vals,occ_check,thresh)
	ctrl_mean = ctrl_diffs.mean(0)
	ctrl_var = ctrl_diffs.var(0)
	for i,j in np.ndindex(np.shape(scan_mean_matrix)):
		exp_times = exp_times_array[i,j]
		exp_vals = exp_vals_array[i,j]
		mean_entry,var_entry = mean_disruption_index(ctrl_mean,ctrl_var,exp_times,exp_vals,occ_check,thresh)
		scan_mean_matrix[i,j] = mean_entry
		scan_var_matrix[i,j] = var_entry
	return scan_mean_matrix,scan_var_matrix

#new
def sequence_deviance_list(ctrl_peak_times,ctrl_peak_vals,exp_times,exp_vals,occ_check=True,thresh=10):
	ctrl_mean,ctrl_var = get_mean_times(ctrl_peak_times,ctrl_peak_vals,thresh,occ_check,False)
	good_exp_times,good_exp_vals = sequence_list_veto(exp_times,exp_vals,occ_check,thresh)
	devs = good_exp_times - ctrl_mean
	devs_std = devs / np.sqrt(ctrl_var)
	#dev_sums = (abs(devs_std)).sum(1) / np.shape(ctrl_mean)[0]
	dev_sums = devs_std.sum(1) / np.shape(ctrl_mean)[0]
	return dev_sums

#new
def scan_deviance(ctrl_peak_times,ctrl_peak_vals,exp_times_array,exp_vals_array,occ_check=True,thresh=10):
	scan_mean_matrix = np.zeros(np.shape(exp_times_array)[0:2])
	scan_var_matrix = np.zeros(np.shape(exp_times_array)[0:2])
	for i,j in np.ndindex(np.shape(scan_mean_matrix)):
		exp_times = exp_times_array[i,j]
		exp_vals = exp_vals_array[i,j]
		devs = sequence_deviance_list(ctrl_peak_times,ctrl_peak_vals,exp_times,exp_vals,occ_check,thresh)
		scan_mean_matrix[i,j] = devs.mean(0)
		scan_var_matrix[i,j] = devs.var(0)
	return scan_mean_matrix,scan_var_matrix

def newclust_sort(w_mat,clust_size=20,clust_n=5):
	n_neurons = np.shape(w_mat)[0]
	first_ext = clust_size * clust_n
	mean_w = w_mat[w_mat > 0].mean()
	mask = w_mat > mean_w
	new_mat_in = np.zeros((n_neurons,n_neurons))
	new_mat_out = np.zeros((n_neurons,n_neurons))
	clust_counts_in = np.zeros((n_neurons,clust_n))
	clust_counts_out = np.zeros((n_neurons,clust_n))
	for i,j in np.ndindex((n_neurons,n_neurons)):
		if w_mat[i,j] > mean_w:
			for k in range(clust_n):
				n_min = k * clust_size
				n_max = (k + 1) * clust_size
				if n_min <= j < n_max:
					clust_counts_out[i,k] += 1
				if n_min <= i < n_max:
					clust_counts_in[j,k] += 1
	max_in_index = np.argmax(clust_counts_in,axis=1)
	max_out_index = np.argmax(clust_counts_out,axis=1)
	new_order_in = np.zeros((n_neurons))
	new_order_out = np.zeros((n_neurons))
	in_counter = 0
	out_counter = 0
	for i in range(clust_n):
		for j in range(n_neurons):
			if (max_in_index[j] == i):
				new_order_in[in_counter] = j
				in_counter += 1
			if (max_out_index[j] == i):
				new_order_out[out_counter] = j
				out_counter += 1
	new_order_in = new_order_in.astype(int)
	new_order_out = new_order_out.astype(int)
	for i,j in np.ndindex((n_neurons,n_neurons)):
		new_mat_in[i,j] = w_mat[new_order_in[i],new_order_in[j]]
		new_mat_out[i,j] = w_mat[new_order_out[i],new_order_out[j]]
	return new_mat_in,new_mat_out




def mean_fb(w_mat,rel=True,el=5,oc=0.5):
#w_mat is recurrent / square matrix
#returns mean feedforward weight
#if rel, relative to mean weight
#assumes 5 element sequence oppcupying first half of network
	n = shape(w_mat)[0]
	n_in = int (n*oc/el)
	w_mean = mean(w_mat[w_mat>0])
	w_ff_list = []
	for i in range(el-1):
		w_ff_temp = w_mat[n_in*(i+1):n_in*(i+2),n_in*(i+1):n_in*(i+2)]
		w_ff_temp_list = w_ff_temp[w_ff_temp>0]
		for w in w_ff_temp_list.T:
			w_ff_list.append(w.item())
	w_ff_mean = mean(w_ff_list)
	if rel:
		return w_ff_mean / w_mean
	else:
		return w_ff_mean

#syns=syn_ee,monitor=weights_ee.w_ampa,num=n_e by default
def mean_fb_array(syns,monitor,num,rel=True):
	#print('a')
	n_t = shape(monitor)[1]
	#print('b')
	mean_fb_list = []
	#print('c')
	for i in range(n_t):
		#print('d1')
		w_temp = array(monsmat(i,syns,monitor,num,sprs=False))
		#print(w_temp)
		#print('d2')
		mean_fb_temp = mean_ff(w_temp,rel=rel)
		#print('d3')
		mean_fb_list.append(mean_fb_temp)
		#print('d4')
	#print('e')
	return array(mean_fb_list)

#new
def binned_spike_count(spike_train,d_t=0.0005*second,start_t=0*second,end_t=350*second,total_t=350*second):
	###MODULARIZE TOTAL TIME
	spike_train = (spike_train / second) * pq.s
	start_t = (start_t / second) * pq.s
	end_t = (end_t / second) * pq.s
	total_t = (total_t / second) * pq.s
	d_t = (d_t / second) * pq.s
	#st = ne.SpikeTrain(spike_train,t_start=start_t,t_stop=end_t)
	#bst = conv.BinnedSpikeTrain(st,binsize=d_t,t_start=start_t,t_stop=end_t)
	st = ne.SpikeTrain(spike_train,t_stop=total_t)
	bst = conv.BinnedSpikeTrain(st,binsize=d_t,t_stop=total_t)
	spike_count = bst.to_array()
	"""
	total_t = end_t - start_t
	total_len = int(total_t / d_t) #+ 1
	spike_count = np.zeros(total_len)
	for i in range(total_len):
		upper = (i + 1) * d_t + start_t
		lower = i * d_t + start_t
		spike_count[i] = len(spike_train[(spike_train > lower) & (spike_train < upper)])
	"""
	return spike_count[0][int(start_t/d_t):int(end_t/d_t)]
	#return spike_count

#new
def binned_spike_counts(spike_monitor,d_t=0.0005*second,start_t=0*second,end_t=350*second):
	#n_neurons = len(spike_monitor.item())
	n_neurons = len(spike_monitor.spike_trains())
	total_t = end_t - start_t
	total_len = int(total_t / d_t) #+ 1
	spike_counts = np.zeros((n_neurons,total_len))
	for i in range(n_neurons):
		#spike_counts[i] = binned_spike_count(spike_monitor.item()[i],d_t=d_t,start_t=start_t,end_t=end_t)
		spike_counts[i] = binned_spike_count(spike_monitor.spike_trains()[i],d_t=d_t,start_t=start_t,end_t=end_t)
	return spike_counts

#new, takes binned_spike_counts with distr. and ctrl. changing over at split_t
def pca_dist_info(spike_counts,d_t=0.0005*second,split_t=100*second, dist_first=True,
		rtrn_spikes=False,traj_window = 0.1*second,cue_freq=2/second,cue_offset=0*second,train='all'):
	#bug fix - WTF
	#cue_freq = 2 / second
	#cue_offset=0*second
	#pdb.set_trace()
	n_trials = int(np.shape(spike_counts)[1] * d_t * cue_freq)
	#print('debug 1:')
	#print(cue_freq)
	#print(d_t)
	cue_steps = int((1 / cue_freq) / d_t)
	window_steps = int(traj_window / d_t)
	offset_steps = int(cue_offset / d_t)
	split_i = int(split_t / d_t)

	if dist_first:
		n_dist_trials = int(split_t * cue_freq)
		n_ctrl_trials = n_trials - n_dist_trials 
	else:
		n_ctrl_trials = int(split_t * cue_freq)
		n_dist_trials = n_trials - n_ctrl_trials


	#ctrl_list = np.zeros(shape=(n_ctrl_trials,np.shape(spike_counts)[0],cue_steps))
	#dist_list = np.zeros(shape=(n_dist_trials,np.shape(spike_counts)[0],cue_steps))
	ctrl_list = np.zeros(shape=(n_ctrl_trials,np.shape(spike_counts)[0],window_steps))
	dist_list = np.zeros(shape=(n_dist_trials,np.shape(spike_counts)[0],window_steps))

	if dist_first:
		for i in range(n_ctrl_trials):
			start_i = (i + n_dist_trials) * cue_steps + offset_steps
			#end_i = (i + n_dist_trials + 1) * cue_steps + offset_steps
			end_i = start_i + window_steps
			temp_st = spike_counts[:,start_i:end_i]
			ctrl_list[i] = temp_st
		for i in range(n_dist_trials):
			start_i = (i) * cue_steps + offset_steps
			#end_i = (i + 1) * cue_steps + offset_steps
			end_i = start_i + window_steps
			temp_st = spike_counts[:,start_i:end_i]
			dist_list[i] = temp_st
	else:
		for i in range(n_dist_trials):
			start_i = (i + n_ctrl_trials) * cue_steps + offset_steps
			#end_i = (i + n_ctrl_trials + 1) * cue_steps + offset_steps
			end_i = start_i + window_steps
			temp_st = spike_counts[:,start_i:end_i]
			dist_list[i] = temp_st
		for i in range(n_ctrl_trials):
			start_i = (i) * cue_steps + offset_steps
			#end_i = (i + 1) * cue_steps + offset_steps
			end_i = start_i + window_steps
			temp_st = spike_counts[:,start_i:end_i]
			ctrl_list[i] = temp_st			

	pca = PCA()
	if train == 'ctrl':
		if dist_first:
			pcmps = pca.fit_transform(spike_counts[split_t:])
			evr = pca.fit(spike_counts[split_t:]).explained_variance_ratio_
		else:
			pcmps = pca.fit_transform(spike_counts[:split_t])
			evr = pca.fit(spike_counts[:split_t]).explained_variance_ratio_
	elif train == 'dis':
		if dist_first:
			pcmps = pca.fit_transform(spike_counts[:split_t])
			evr = pca.fit(spike_counts[:split_t]).explained_variance_ratio_
		else:
			pcmps = pca.fit_transform(spike_counts[split_t:])
			evr = pca.fit(spike_counts[split_t:]).explained_variance_ratio_
	else:
		pcmps = pca.fit_transform(spike_counts)
		evr = pca.fit(spike_counts).explained_variance_ratio_

	traj_ctrl_list = np.zeros(shape=np.shape(ctrl_list))
	for i in range(len(ctrl_list)):
		traj_ctrl_list[i] = np.dot(pcmps,ctrl_list[i])
	traj_dist_list = np.zeros(shape=np.shape(dist_list))
	for i in range(len(dist_list)):
		traj_dist_list[i] = np.dot(pcmps,dist_list[i])		

	return pcmps,evr,traj_ctrl_list,traj_dist_list

def pca_traj_diff(traj_ctrl_list,traj_dist_list):
	mean_traj_ctrl = traj_ctrl_list.mean(0)
	#std_traj_ctrl = traj_ctrl_list.std(0)
	traj_diffs = np.zeros(shape=np.shape(traj_dist_list))
	for i in range(np.shape(traj_dist_list)[0]):
		traj_diffs[i] = traj_dist_list[i] - mean_traj_ctrl
	mean_traj_diffs = traj_diffs.mean(0)
	std_traj_diffs = traj_diffs.std(0)
	return mean_traj_diffs,std_traj_diffs

def pca_traj_dist(evr,traj_ctrl_list,traj_dist_list): #"evr normed euclidean distance"
	mean_traj_diffs,std_traj_diffs = pca_traj_diff(traj_ctrl_list,traj_dist_list)
	traj_dist = np.zeros(shape=np.shape(traj_ctrl_list)[2])
	for i in range(np.shape(evr)[0]):
		traj_dist += (evr[i] * mean_traj_diffs[i])**2
	return np.sqrt(traj_dist)

def pca_traj_dist_noweight(traj_ctrl_list,traj_dist_list): #"evr normed euclidean distance"
	mean_traj_diffs,std_traj_diffs = pca_traj_diff(traj_ctrl_list,traj_dist_list)
	traj_dist = np.zeros(shape=np.shape(traj_ctrl_list)[2])
	for i in range(np.shape(traj_dist_list)[1]):
		traj_dist += (mean_traj_diffs[i])**2
	return np.sqrt(traj_dist)

def ptd_from_st(spike_monitor,d_t=0.0005*second,start_t=150*second,end_t=350*second,split_t=100*second,
		dist_first=True,traj_window = 0.1*second,cue_freq=2/second,cue_offset=0*second,train='all'):
	spike_counts = binned_spike_counts(spike_monitor,d_t,start_t,end_t)
	#print('debug 2:')
	#print(cue_freq)
	#print(d_t)
	pc,ev,tc,td = pca_dist_info(spike_counts,d_t,split_t,dist_first,False,traj_window,cue_freq,cue_offset,train)
	ptd = pca_traj_dist(ev,tc,td)
	return ptd

def ptd_from_st_noweight(spike_monitor,d_t=0.0005*second,start_t=150*second,end_t=350*second,split_t=100*second,
		dist_first=True,traj_window = 0.1*second,cue_freq=2/second,cue_offset=0*second,train='all'):
	spike_counts = binned_spike_counts(spike_monitor,d_t,start_t,end_t)
	#print('debug 2:')
	#print(cue_freq)
	#print(d_t)
	pc,ev,tc,td = pca_dist_info(spike_counts,d_t,split_t,dist_first,False,traj_window,cue_freq,cue_offset,train)
	#pdb.set_trace()
	ptd = pca_traj_dist_noweight(tc,td)
	return ptd

def ptd_exp_array(exp_array,d_t=0.0005*second,start_t=150*second,end_t=350*second,split_t=100*second,
		dist_first=True,traj_window = 0.1*second,cue_freq=2/second,cue_offset=0*second,train='all'):
	with np.errstate(divide='ignore', invalid='ignore'): #cheat my way through this
		exp_shape = np.shape(exp_array)
		window_steps = int(traj_window / d_t)
		#ptd_array = zeros(shape=(4,4,5,window_steps)) #CHEATING - AUTOMATE LATER
		ptd_array = zeros(shape=(exp_shape[0],exp_shape[1],exp_shape[2],window_steps))
		for i,j,k in ndindex(exp_shape):
			print('exp_array element '+str(i)+', '+str(j)+', '+str(k))
			#print('debug 3:')
			#print(cue_freq)
			#print(d_t)
			ptd_array[i][j][k] = ptd_from_st(exp_array[i][j][k],d_t,start_t,end_t,split_t,dist_first,traj_window,cue_freq,cue_offset,train)
			gc.collect()
	return ptd_array

def ptd_exp_array_noweight(exp_array,d_t=0.0005*second,start_t=150*second,end_t=350*second,split_t=100*second,
		dist_first=True,traj_window = 0.1*second,cue_freq=2/second,cue_offset=0*second,train='all'):
	with np.errstate(divide='ignore', invalid='ignore'): #cheat my way through this
		exp_shape = np.shape(exp_array)
		window_steps = int(traj_window / d_t)
		#ptd_array = zeros(shape=(4,4,5,window_steps)) #CHEATING - AUTOMATE LATER
		ptd_array = zeros(shape=(exp_shape[0],exp_shape[1],exp_shape[2],window_steps))
		for i,j,k in ndindex(exp_shape):
			print('exp_array element '+str(i)+', '+str(j)+', '+str(k))
			#print('debug 3:')
			#print(cue_freq)
			#print(d_t)
			ptd_array[i][j][k] = ptd_from_st_noweight(exp_array[i][j][k],d_t,start_t,end_t,split_t,dist_first,traj_window,cue_freq,cue_offset,train)
			gc.collect()
	return ptd_array


#HIGHLY IMPROVED
def simple_decision_readoutv2(rate_list,min_time=150.,max_time=250.,cue_freq=2.,window=0.035,offset=-0.01):
	dt = 10000
	active = 4
	trials = int(round((max_time - min_time) * cue_freq))
	starts = np.arange(min_time + offset,max_time + offset,(1 / cue_freq))
	stops = np.arange(min_time + offset + window,max_time + offset + window,(1 / cue_freq))
	peak_time_list = zeros((trials,active))
	peak_val_list = zeros((trials,active)) 
	t_count = 0
	for cue in range(trials):
		start = starts[cue]
		stop = stops[cue]
		traces = rate_list[:,int(round(start*dt)):int(round(stop*dt))]
		for cluster in range(active):
			trace = traces[cluster]
			peak_pos = get_peaks(trace,thresh=0,subsamp=1,min_dist=0)
			if np.size(peak_pos) < 1:
				peak_pos = 999 #failure flag
				peak_val = 0 #failute flag
			elif np.size(peak_pos) > 1:
				mult_peak_vals = trace[np.array(np.array(peak_pos) * dt,dtype=int)]
				peak_index = np.argmax(mult_peak_vals)
				peak_pos = peak_pos[peak_index]
				peak_val = trace[int(round(peak_pos * dt))]
			else:
				peak_pos = peak_pos[0]
				peak_val = trace[int(round(peak_pos * dt))]
			peak_time_list[t_count,cluster] = peak_pos
			peak_val_list[t_count,cluster] = peak_val
		t_count += 1
	decisions = []
	decision_strengths = zeros((trials))
	for trial in range(trials):
		decision = 'err'
		if peak_time_list[trial,0] > peak_time_list[trial,1]:
			decision = 'nb'
		elif (peak_time_list[trial,2] < peak_time_list[trial,1]) and (peak_time_list[trial,3] < peak_time_list[trial,1]):
			decision = 'ncd'
		else:
			if peak_val_list[trial,2] > peak_val_list[trial,3]:
				decision = 'c'
			elif peak_val_list[trial,3] > peak_val_list[trial,2]:
				decision = 'd'
			decision_strengths[trial] = abs(peak_val_list[trial,2] - peak_val_list[trial,3])
		decisions.append(decision)

	return peak_time_list,peak_val_list,decision_strengths,decisions

def get_decision_stacksv2(rate_lists,min_time=150.,max_time=250.,cue_freq=2.,window=0.035,offset=-0.01):
	numlists = shape(rate_lists)[0]
	ds_total = []
	ds_stack = []
	dc_total = []
	dc_stack = []
	for i in range(numlists):
		pt,pv,ds,dc = simple_decision_readoutv2(rate_lists[i,0:4,:],min_time,max_time,cue_freq,window,offset)
		ds_stack.append(ds)
		dc_stack.append(dc)
	for i in range(numlists):
		for j in range(len(dc_stack[i])):
			ds_total.append(ds_stack[i][j])
			dc_total.append(dc_stack[i][j])
	return dc_total,dc_stack,ds_total,ds_stack

def get_good_strengths(dss,dcs):
	numruns = shape(dss)[0]
	dsl = []
	for i in range(numruns):
		for j in range(len(dss[i])):
			if (dcs[i][j] == 'c') or (dcs[i][j] == 'd'):
				dsl.append(dss[i][j])
	return dsl

def rep_weight_sort(weights,clust_size=20,active_clust_n=5):
	all_weights = []
	rec_weights = []
	one_forward = []
	n_forward = []
	one_backward = []
	n_backward = []
	to_ext = []
	from_ext = []
	n_neurons = np.shape(weights)[0]
	total_clusts = int(n_neurons / clust_size)
	for ci in range(total_clusts):
		for cj in range(total_clusts):
			for ni in range(clust_size):
				for nj in range(clust_size):
					i = ci * clust_size + ni
					j = cj * clust_size + nj
					if weights[i,j] > 0:
						current_weight = weights[i,j]
						all_weights.append(current_weight)
						if (ci < (active_clust_n - 1)) and (cj == (ci + 1)):
							one_forward.append(current_weight)
						if (ci < active_clust_n) and (ci == cj):
							rec_weights.append(current_weight)
						if (ci < (active_clust_n - 2)) and (cj > (ci+1)) and (cj < active_clust_n):
							n_forward.append(current_weight)
						if (ci < active_clust_n) and (cj == (ci - 1)):
							one_backward.append(current_weight)
						if (ci < active_clust_n) and (cj < (ci-1)):
							n_backward.append(current_weight)
						if (ci < active_clust_n) and (cj >= active_clust_n):
							to_ext.append(current_weight)
						if (ci >= active_clust_n) and (cj < active_clust_n):
							from_ext.append(current_weight)
	return all_weights,rec_weights,one_forward,n_forward,one_backward,n_backward,to_ext,from_ext

def rep_weights_from_stack(weight_stack,clust_size=20,active_clust_n=5):
	all_weights = []
	rec_weights = []
	one_forward = []
	n_forward = []
	one_backward = []
	n_backward = []
	to_ext = []
	from_ext = []
	for i in range(len(weight_stack)):
		taw,trw,tof,tnf,tob,tnb,tte,tfe = rep_weight_sort(weight_stack[i],clust_size,active_clust_n)
		for w in taw:
			all_weights.append(w)
		for w in trw:
			rec_weights.append(w)
		for w in tof:
			one_forward.append(w)
		for w in tnf:
			n_forward.append(w)
		for w in tob:
			one_backward.append(w)
		for w in tnb:
			n_backward.append(w)
		for w in tte:
			to_ext.append(w)
		for w in tfe:
			from_ext.append(w)
	return all_weights,rec_weights,one_forward,n_forward,one_backward,n_backward,to_ext,from_ext