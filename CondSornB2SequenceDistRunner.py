from CondSornB2ParSequence import *
#import SequenceAnalyzers as SA
import pdb
import os
import gc

def run_it(tag,ei,ie,dist_num,dist_delay,dist_frac=1,trig_frac=1,ctrl=False,play=False,save_it=True):
    dense_ei = ei
    dense_ie = ie
    ### Neuron model
    #
    # We use a noisy conductance based model with excitatory / ampa and inhibitory /
    #  gaba conductances as well as additional variables to permit implementation of
    #  synaptic normalization. We use two versions, one for e and one for i.
    #
    eqs_neuron_e='''
    dv/dt=(g_leak*(v_rest-v)+i_ext+i_syn)/c+s_n*xi*tau_n**-0.5: volt  # voltage
    i_ext : amp  # external current
    i_syn=g_ampa*(e_ampa-v)+g_gaba*(e_gaba-v) : amp  # synaptic current
    dg_ampa/dt=-g_ampa/tau_ampa : siemens  # ampa synaptic conductance
    sumw_ampa : siemens  # total ampa input
    sumw_ampa_target : siemens  # target total ampa input
    dg_gaba/dt=-g_gaba/tau_gaba : siemens  # gaba synaptic conductance
    sumw_gaba : siemens  # toal gaba conductance
    sumw_gaba_target : siemens  # target total gaba input
    dv_t/dt=-eta_ip_decay : volt # firing threshold with linear decay
    refrac : second  # refractory period
    '''
    eqs_neuron_i='''
    dv/dt=(g_leak*(v_rest-v)+i_ext+i_syn)/c+s_n*xi*tau_n**-0.5: volt  # voltage
    i_ext : amp  # external current
    i_syn=g_ampa*(e_ampa-v)+g_gaba*(e_gaba-v) : amp  # synaptic current
    dg_ampa/dt=-g_ampa/tau_ampa : siemens  # ampa synaptic conductance
    sumw_ampa : siemens  # total ampa input
    sumw_ampa_target : siemens  # target total ampa input
    dg_gaba/dt=-g_gaba/tau_gaba : siemens  # gaba synaptic conductance
    sumw_gaba : siemens  # toal gaba conductance
    sumw_gaba_target : siemens  # target total gaba input
    dv_t/dt=-eta_ip_decay : volt # firing threshold with linear decay
    refrac : second  # refractory period
    '''

    ### Synapse model
    #
    # We use conductance based STDP ampa synapses and conductance based gaba synapses
    #  with additional variables for synaptic normalization (not yet implemented). 
    #
    #### DOUBLE CHECK THAT YOU HAVE SUMMED PRE / POST RIGHT
    ampa_model_eq='''
    w_ampa : siemens  # synaptic weight (ampa synapse)
    dwpre/dt=-wpre/taupre : siemens (event-driven)
    dwpost/dt=-wpost/taupost : siemens (event-driven)
    sumw_ampa_post = w_ampa : siemens (summed)
    plastic: 1
    '''
    simple_ampa_model_eq='''
    w_ampa : siemens  # synaptic weight (ampa synapse)
    '''
    ampa_pre_eq='''
    g_ampa_post+=w_ampa  # increment ampa conductance
    wpre+=(dwpre*plastic)
    w_ampa=clip(w_ampa+(wpost*plastic),0,w_max)  # hard boundaries
    '''
    simple_ampa_pre_eq='''
    g_ampa_post+=w_ampa  # increment ampa conductance
    '''
    ampa_post_eq='''
    wpost+=dwpost
    w_ampa*=(sumw_ampa_target/sumw_ampa)
    w_ampa=clip(w_ampa+wpre,0,w_max)  # hard boundaries
    '''
    gaba_model_eq='''
    w_gaba : siemens  # synaptic weight (gaba synapse)
    sumw_gaba_post = w_gaba : siemens (summed)
    '''
    gaba_pre_eq='''
    g_gaba_post+=w_gaba  # increment gaba conductance
    '''
    ### Network model
    #
    # We base the network design on the network described in Miner and Triesch PLOS
    #  CB 2016, minus the spatial topology.
    #
    reset_eqs_e='''
    v=v_rest  # reset to resting potential upon spike
    v_t+=eta_ip_spike  # increase threshold voltage
    '''
    reset_eqs_i='''
    v=v_rest  # reset to resting potential upon spike
    v_t+=(eta_ip_spike/5.)  # increase threshold voltage
    '''
    """
    ##commented out to account for subgroup bug with summed variables, will return
    # to single group w/ subgroups when fixed; further corrections propagated through
    # code
    #neurons=NeuronGroup(n_e+n_i,eqs_neuron,threshold=thr_eqs,reset=reset_eqs)
    neurons=NeuronGroup(n_e+n_i,eqs_neuron,threshold='v>v_t',reset=reset_eqs,
                        refractory='refrac')
    nrns_e=neurons[:n_e]  # excitatory subgroup
    nrns_i=neurons[n_e:]  # inhibitory subgroup
    """
    nrns_e=NeuronGroup(n_e,eqs_neuron_e,threshold='v>v_t',reset=reset_eqs_e,
                        refractory='refrac')
    nrns_i=NeuronGroup(n_i,eqs_neuron_i,threshold='v>v_t',reset=reset_eqs_i,
                        refractory='refrac')
    syn_ee = Synapses(nrns_e,nrns_e,model=ampa_model_eq,
                    on_pre=ampa_pre_eq,on_post=ampa_post_eq,delay=d_ee)  # create e->e synapses
    #syn_ee = Synapses(nrns_e,nrns_e,model=ampa_model_eq,
    #                on_pre=simple_ampa_pre_eq)  # create non-plastic e->e synapses
    syn_ee.connect(condition='i!=j',p=dense_ee)  # connect e->e synapses
    syn_ei = Synapses(nrns_e,nrns_i,model=ampa_model_eq,
                    on_pre=simple_ampa_pre_eq,delay=d_ei)  # create e->i synapses
    syn_ei.connect(p=dense_ei)  # connect e->i synapses
    syn_ie = Synapses(nrns_i,nrns_e,model=gaba_model_eq,
                    on_pre=gaba_pre_eq,delay=d_ie)  # create i->e synapses
    syn_ie.connect(p=dense_ie)  # connect i->e synapses
    syn_ii = Synapses(nrns_i,nrns_i,model=gaba_model_eq,
                    on_pre=gaba_pre_eq,delay=d_ii) # create i->i synapses
    syn_ii.connect(condition='i!=j',p=dense_ii)  # connect i->i synapses
    ##init hack
    #w_start=0.15*nS
    w_start=0.5*nS
    syn_ee.w_ampa=w_start
    syn_ei.w_ampa=w_start*2.0
    syn_ie.w_gaba=w_start*2.0
    syn_ii.w_gaba=w_start
    nrns_e.sumw_ampa_target = w_start * n_e * dense_ee
    #neurons.v=np.random.normal(2.5+v_rest/mV,2.5,n_all)*mV
    nrns_e.v=np.random.normal(2.5+v_rest/mV,2.5,n_e)*mV
    nrns_i.v=np.random.normal(2.5+v_rest/mV,2.5,n_i)*mV
    #neurons.v=v_rest
    #neurons.v_t=np.random.normal(2.5+v_thresh/mV,5,n_all)*mV
    #neurons.v_t=np.random.normal(7.5+v_rest/mV,2.5,n_all)*mV
    nrns_e.v_t=np.random.normal(3.5+v_rest/mV,2.5,n_e)*mV
    nrns_i.v_t=np.random.normal(3.5+v_rest/mV,2.5,n_i)*mV
    #neurons.v_t=v_thresh
    nrns_e.refrac=refrac_e
    nrns_i.refrac=refrac_i
    ##ADD DELAYS later

    ### Stimulation
    stim_rate = 50
    cluster_size = 20
    stim = TimedArray(tile([[stim_rate,0.,0.,0.,0.,0.,0.,0.,0.,0.],
        [0.,stim_rate,0.,0.,0.,0.,0.,0.,0.,0.],
        [0.,0.,stim_rate,0.,0.,0.,0.,0.,0.,0.],
        [0.,0.,0.,stim_rate,0.,0.,0.,0.,0.,0.],
        [0.,0.,0.,0.,stim_rate,0.,0.,0.,0.,0.]]*Hz,
        int((total_time)/second)).T,dt=0.1*second)
    # set up input subgroups
    # random at the moment (except by index); spatial clustering later?
    gin1 = nrns_e[0:cluster_size]
    gin2 = nrns_e[cluster_size:2*cluster_size]
    gin3 = nrns_e[2*cluster_size:3*cluster_size]
    gin4 = nrns_e[3*cluster_size:4*cluster_size]
    gin5 = nrns_e[4*cluster_size:5*cluster_size]
    ###SPECIFIC TO POP 200 config, remove later
    gin6 = nrns_e[5*cluster_size:6*cluster_size]
    gin7 = nrns_e[6*cluster_size:7*cluster_size]
    gin8 = nrns_e[7*cluster_size:8*cluster_size]
    gin9 = nrns_e[8*cluster_size:9*cluster_size]
    gin10 = nrns_e[9*cluster_size:10*cluster_size]
    ###SPECIFIC TO POP 200 config, remove later
    # set up poisson inputs
    pin = PoissonGroup(5, rates='stim(t,i)')
    #pin = PoissonGroup(5, rates=lambda t: stim(t))
    #MSpin = SpikeMonitor(pin)
    # set up connection to main network
    syn_in = Synapses(pin,nrns_e,model=simple_ampa_model_eq,on_pre=simple_ampa_pre_eq) 
    #cin = Connection(pin, nrns_e, 'g_ampa')
    #syn_in.connect(i=pin[0:1],j=gin1)
    #syn_in.connect(i=pin[1:2],j=gin2)
    #syn_in.connect(i=pin[2:3],j=gin3)
    #syn_in.connect(i=pin[3:4],j=gin4)
    #syn_in.connect(i=pin[4:5],j=gin5)
    syn_in.connect(i=0,j=arange(0,cluster_size))
    syn_in.connect(i=1,j=arange(cluster_size,2*cluster_size))
    syn_in.connect(i=2,j=arange(2*cluster_size,3*cluster_size))
    syn_in.connect(i=3,j=arange(3*cluster_size,4*cluster_size))
    syn_in.connect(i=4,j=arange(4*cluster_size,5*cluster_size))
    syn_in.w_ampa = 0.0*nS
    #syn_in.w_ampa = 0.25*nS
    #gtrigger = nrns_e[0:cluster_size]
    #syn_trig = Synapses(pin,nrns_e,model=ampa_model_eq,on_pre=simple_ampa_pre_eq) 
    #syn_trig.connect(i=0, j=arange(0,cluster_size))
    #syn_trig.w_ampa = 0.0*nS

    pin_trig = SpikeGeneratorGroup(1,np.array([0]),np.array([0])*second,period=0.5*second)
    #gtrigger = nrns_e[0:(int(cluster_size/4))]
    syn_trig = Synapses(pin_trig,nrns_e,model=simple_ampa_model_eq,on_pre=simple_ampa_pre_eq)
    trig_size = int(cluster_size/trig_frac)
    syn_trig.connect(i=0, j=arange(0,trig_size))
    syn_trig.w_ampa = 0.0*nS

    #dist_delay = 0.000 #in seconds, total replay in default params takes ~0.006
    pin_dist = SpikeGeneratorGroup(1,np.array([0]),np.array([0+dist_delay])*second,period=0.5*second)
    #gdist = nrns_e[0:(int(cluster_size/4))]
    syn_dist = Synapses(pin_dist,nrns_e,model=simple_ampa_model_eq,on_pre=simple_ampa_pre_eq)
    #syn_dist.connect(i=0, j=arange(0,cluster_size/2))
    #dist_num = 2
    #dist_frac = 2
    dist_size = int(cluster_size/dist_frac)
    syn_dist.connect(i=0, j=arange((dist_num-1)*cluster_size,(((dist_num-1)*cluster_size)+dist_size)))
    syn_dist.w_ampa = 0.0*nS

    ### Monitors
    defaultclock.dt=0.1*ms

    #spikes=SpikeMonitor(neurons)
    spikes_e=SpikeMonitor(nrns_e)
    spikes_i=SpikeMonitor(nrns_i)

    neuron_mon_var=['v','g_ampa','g_gaba','v_t','sumw_ampa']
    #neuron_mon=StateMonitor(neurons,variables=neuron_mon_var,record=True,dt=10.0*ms)
    neuron_mon_e=StateMonitor(nrns_e,variables=neuron_mon_var,record=True,dt=10.0*ms)
    neuron_mon_i=StateMonitor(nrns_i,variables=neuron_mon_var,record=True,dt=10.0*ms)
    rate_e=PopulationRateMonitor(nrns_e)
    rate_i=PopulationRateMonitor(nrns_i)
    MRgin1=PopulationRateMonitor(gin1)
    MRgin2=PopulationRateMonitor(gin2)
    MRgin3=PopulationRateMonitor(gin3)
    MRgin4=PopulationRateMonitor(gin4)
    MRgin5=PopulationRateMonitor(gin5)
    MRgin6=PopulationRateMonitor(gin6)
    MRgin7=PopulationRateMonitor(gin7)
    MRgin8=PopulationRateMonitor(gin8)
    MRgin9=PopulationRateMonitor(gin9)
    MRgin10=PopulationRateMonitor(gin10)
    rate_list=[MRgin1,MRgin2,MRgin3,MRgin4,MRgin5,MRgin6,MRgin7,MRgin8,MRgin9,MRgin10]

    weights_ee=StateMonitor(syn_ee,'w_ampa',record=True,dt=100.0*ms)
    net=Network(collect())
    #pdb.set_trace() #XXXX
    print("Warm Up")
    syn_ee.plastic=1
    net.run(warm_time,report='text')
    net.store('warmed_up')
    net.restore('warmed_up')
    print("Training")
    syn_in.w_ampa = 20.0*nS
    net.run(train_time,report='text')
    net.store('trained')
    net.restore('trained')
    print("Relaxing")
    syn_in.w_ampa = 0.0*nS
    syn_ee.plastic=0
    net.run(relax_time,report='text')
    net.store('relaxed')
    net.restore('relaxed')
    print("Testing Phase 1")
    syn_trig.w_ampa = 20.0*nS
    if not ctrl:
        syn_dist.w_ampa = 20.0*nS
    net.run(test_time,report='text')
    net.store('test1')
    net.restore('test1')
    print("Testing Phase 2")
    syn_dist.w_ampa = 0.0*nS
    net.run(test_time,report='text')


    #helper
    def monasmat(k,syns=syn_ee,monitor=weights_ee.w_ampa,num=n_e,sprs=True):
        temp_mat = zeros((num,num))
        for i,j,w in zip(syns.i,syns.j,monitor[:,k]):
            temp_mat[i,j] = w / siemens
        return temp_mat 
        """
        if sprs:
            temp_mat = sparse.lil_matrix((num,num))
            for i,j,w in zip(syns.i,syns.j,monitor[:,k]):
                temp_mat[i,j] = w / siemens
            return temp_mat
        else:
            temp_mat = zeros((num,num))
            for i,j,w in zip(syns.i,syns.j,monitor[:,k]):
                temp_mat[i,j] = w / siemens
            return temp_mat
        """

    if play:
        interactive = True
        #anim_w = True
    else:
        interactive = False
    anim_w = False
    if interactive:
        figure()
        subplot(211)
        plot(rate_e.t,rate_e.smooth_rate(window='gaussian',width=250.*ms)/Hz)
        subplot(212)
        plot(rate_i.t,rate_i.smooth_rate(window='gaussian',width=250.*ms)/Hz)
        figure()
        plot(MRgin1.t,MRgin1.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="A")
        plot(MRgin2.t,MRgin2.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="B")
        plot(MRgin3.t,MRgin3.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="C")
        plot(MRgin4.t,MRgin4.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="D")
        plot(MRgin5.t,MRgin5.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="E")
        plot(MRgin10.t,MRgin6.smooth_rate(window='gaussian',width=2.*ms)/Hz,label="CONTROL")
        ratemon_list = [MRgin1,MRgin2,MRgin3,MRgin4,MRgin5,MRgin6,MRgin10]
        legend()
        figure()
        subplot(311)
        #plot(spikes.t,spikes.i,'.')
        plot(spikes_e.t,spikes_e.i,'.')
        plot(spikes_i.t,spikes_i.i,'.')
        subplot(312)
        plot(neuron_mon_e.t,neuron_mon_e.v[0],'-')
        plot(neuron_mon_e.t,neuron_mon_e.v[1],'-')
        plot(neuron_mon_e.t,neuron_mon_e.v[2],'-')
        plot(neuron_mon_e.t,neuron_mon_e.v[3],'-')
        plot(neuron_mon_e.t,neuron_mon_e.v[4],'-')
        plot(neuron_mon_e.t,neuron_mon_e.v_t[0],'--')
        plot(neuron_mon_e.t,neuron_mon_e.v_t[1],'--')
        plot(neuron_mon_e.t,neuron_mon_e.v_t[2],'--')
        plot(neuron_mon_e.t,neuron_mon_e.v_t[3],'--')
        plot(neuron_mon_e.t,neuron_mon_e.v_t[4],'--')
        subplot(313)
        """
        plot(neuron_mon.t,neuron_mon.v[200],'-')
        plot(neuron_mon.t,neuron_mon.v[201],'-')
        plot(neuron_mon.t,neuron_mon.v[202],'-')
        plot(neuron_mon.t,neuron_mon.v[203],'-')
        plot(neuron_mon.t,neuron_mon.v[204],'-')
        plot(neuron_mon.t,neuron_mon.v_t[200],'--')
        plot(neuron_mon.t,neuron_mon.v_t[201],'--')
        plot(neuron_mon.t,neuron_mon.v_t[202],'--')
        plot(neuron_mon.t,neuron_mon.v_t[203],'--')
        plot(neuron_mon.t,neuron_mon.v_t[204],'--')
        """
        plot(neuron_mon_i.t,neuron_mon_i.v[0],'-')
        plot(neuron_mon_i.t,neuron_mon_i.v[1],'-')
        plot(neuron_mon_i.t,neuron_mon_i.v[2],'-')
        plot(neuron_mon_i.t,neuron_mon_i.v[3],'-')
        plot(neuron_mon_i.t,neuron_mon_i.v[4],'-')
        plot(neuron_mon_i.t,neuron_mon_i.v_t[0],'--')
        plot(neuron_mon_i.t,neuron_mon_i.v_t[1],'--')
        plot(neuron_mon_i.t,neuron_mon_i.v_t[2],'--')
        plot(neuron_mon_i.t,neuron_mon_i.v_t[3],'--')
        plot(neuron_mon_i.t,neuron_mon_i.v_t[4],'--')
        show()
        pdb.set_trace()
    if anim_w:
        import matplotlib.animation as animation
        fig = figure()
        ims = []
        #for i in range(int((total_time)/(100.0*ms))):
            #im = imshow(monasmat(k=i,sprs=False),cmap='cubehelix')
            #colorbar()
            #ims.append([im])
        #ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,repeat_delay=1000)
        for i in range(int((total_time)/(100.0*ms))):
            ims.append([imshow(monasmat(k=i,sprs=False),cmap='cubehelix')]) #change time
        ani = animation.ArtistAnimation(fig, ims, interval=5, blit=True,repeat_delay=1000)
        ani.save('animation3.mp4')
        show()
        pdb.set_trace()
    if save_it:
        savepath = "data"
        fname = "w_ee_final_"+tag+".npy"
        w_mat = monasmat(k=int((total_time/second)*10-1),sprs=False)
        save(os.path.join(savepath,fname),w_mat)
        fname = "spike_trains_"+tag+".npz"
        spiketrains = spikes_e.spike_trains()
        save(os.path.join(savepath,fname),spiketrains)
        ratemon_list = [MRgin1,MRgin2,MRgin3,MRgin4,MRgin5,MRgin6,MRgin10]
        rate_array_list = []
        #2ms window; find more flexible solution?
        for r_m in ratemon_list:
            rate_array_list.append(r_m.smooth_rate(window='gaussian',width=2.*ms)/Hz)
        fname = "rate_list"+tag+".npy"
        save(os.path.join(savepath,fname),rate_array_list)
        #save(os.path.join(savepath,fname),ratemon_list)
    gc.collect()