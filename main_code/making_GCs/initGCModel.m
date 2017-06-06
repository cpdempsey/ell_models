function GC_model = initGCModel

%this gets called when we launch the gui or the batch fitter, so we can
%load all the data and initialize the model. if you want to make changes to
%the default model settings, make them here!

GC_model=struct;
GC_model.tau_m       = 8.7; % 8.7 (ms)
GC_model.v_thresh    = -43; % (mV);
GC_model.rmgs        = 1; % unitless
GC_model.El         = -63; % (mV)
GC_model.V_reset     = GC_model.El; % (mV);
GC_model.tau_s       = .2; % (ms) fast synaptic conductance
GC_model.tau_s_slow  = 37.8; %(ms) slow synaptic conductance
GC_model.Ws          = GC_model.rmgs*zeros(1,6); 
GC_model.spiking_on  = 0;

GC_model.mf_input    = [0 0 0];
GC_model.MF_prob     = [1 1 1]; %probability of MF spiking on a trial
GC_model.GC_to_model = 1;

GC_model.dt          = 5e-2; %Nate uses 5e-2; we need to step this up to 5e-3.
GC_model.min_t       = -.025e3;
GC_model.max_t       = .2e3;

GC_model.tRefrac     = 4/GC_model.dt; %add a 4ms refractory period
GC_model.EPSPvar     = 0.05; %this is the average EPSP variance in the cells we looked at.

GC_model             = calc_EPSP_scale(GC_model); %gets the amplitude scaling of a spike by the synapse

end