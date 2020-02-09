from brian2 import *

### Base parameters
warm_time = 50 * second
train_time = 50 * second
relax_time = 50 * second
test_time = 100 * second
total_time = warm_time + train_time + relax_time + test_time # run time
#slice_x = 1000.0  #virtual slice x, microns
#slice_y = 1000.0  #virtual slice y, microns
#con_type = 'norm' #connection profile type
#con_type = 'unif' #uniform connection profile; widths ignored
##width_ei = 200  #ei connectivity radius, microns
#width_ie = 200  #ie connectivity radius, microns
#width_ii = 200  #ii connectivity radius, microns
dense_ee = 0.2 # recurrent excitatory sparseness
dense_ie = 0.2 # inhibitory to excitatory sparseness
dense_ei = 0.2 #  excitatory to inhibitory  sparseness
dense_ii = 0.0 # inhibitory to inhibitory sparseness
#dense_ii = 0.5# inhibitory to inhibitory sparseness
n_e=200
n_i=int(0.2*n_e)
n_all=n_e+n_i
"""
d_ee = 1.5 * ms # e->e latency for 1 mm
d_ei = 0.5 * ms # e->i latency for 1 mm
d_ie = 1.0 * ms # i->e latency for 1 mm
d_ii = 1.0 * ms # i->i latency for 1 mm
"""
d_ee = 0. * ms # e->e latency for 1 mm
d_ei = 0. * ms # e->i latency for 1 mm
d_ie = 0. * ms # i->e latency for 1 mm
d_ii = 0. * ms # i->i latency for 1 mm

### Neuron and spike plasticity parameters
g_leak = 30. * nS # leak conductance
v_rest = -70. * mV # resting potential
c = 300. * pF # membrane capacitance
s_n = sqrt(1.0) * mV # noise amplitude
tau_n = 20. * ms # noise timescale
tau_lp1 = 40. * ms # first voltage lp filter timescale
tau_lp2 = 30. * ms # second voltage lp filter timescale
tau_lph = 1000. * ms # homeostatic metaplasticity voltage lp filter timescale
e_ampa = 0. * mV # ampa reversal potential
e_gaba = -85. * mV # gaba reversal potential
tau_ampa = 2. * ms # ampa timescale
tau_gaba = 5. * ms # gaba timescale
tau_y = 33. * ms # homeostatic intrinsic plasticity spike trace timsecale
y_trgt_df = 0.5 # default target y value
eta_vt = 1. * mV / second # homeostatic intrinsic plasticity adaptation rate
refrac_e = 10. * ms # excitatory refractory period
refrac_i = 2. * ms # inhibitory refractory period
#refrac_e = 2. * ms 
#refrac_i = 1. * ms 
####
taupre = 20 * ms # XXX
taupost = taupre # XXX
dwpre = 5e-2 * nS # XXX
dwpost = -dwpre*taupre/taupost # XXX
w_max = 50. * nS # maximum ampa conductance weight

#v_thresh = -55. * mV # firing threshold
v_thresh = -66. * mV
y_reset = .1 # ip spike trace reset value

### IP parameters
eta_ip_decay = 0.2 * mV / second # IP trace decay rate
eta_ip_spike = 0.066 * mV # IP addition on spike

### Synaptic normalization parameters
#c_scale = 4.0
#c_scale = 5.0
#norm_rate_ee = 1.0
#norm_ee = n_e * dense_ee * w_ee * c_scale
#norm_ie = n_i * dense_ie * w_ie * c_scale #* 2.0
#norm_ei = n_e * dense_ei * w_ei * c_scale #* 2.0
#norm_ii = n_i * dense_ii * w_ii * c_scale

### Stimulation parameters
#stim_rate = 100 #stim rate in Hz
#stim_w = 1.0 #stim connection weight in unitless conductance
#trig_w = 0.25
#cluster_size = 20 #cluster size
#num_clusts = 5 #number of clusters, just a specifier - not automatic (yet)

### Finishing elements
#file_tag = 'scratch'
anim_w = False
#save_w_final = False
#save_w_full = False
interactive = True
