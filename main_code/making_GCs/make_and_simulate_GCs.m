cd('~/Documents/MATLAB/ELL_network_project/');

clear all

% where to save generated GC data
savepth = './';
expdir  = 'generated_GCs';
expname = expdir;
if(~exist([savepth expdir],'dir'))
    mkdir(savepth,expdir);
end

% file containing mossy fiber data used as inputs to GCs
mossy_file = './data/mossy_resps.mat';

% load the mossy fibers
mfinp           =   load(mossy_file,'mftypes','rspstore');
rspstore        =   mfinp.rspstore; % times
mftypes         =   mfinp.mftypes;  % classification of MF


% where are the fit GC parameters and indices to fit mossy fibers
GCparams_file = './data/GCparams_file.mat';

% load Ann's parameters:
loadGCP=load(GCparams_file,'MF_indices','GCparams');
GCparams=loadGCP.GCparams;
MF_indices=loadGCP.MF_indices;

GCparams.Wsparse(170,:)=0; %that last pause GC is trouble!

%% make the granule cells (this will take a minute)

N = 20000;
[GC_models, celltypes] = generateGCPopAnn(N,rspstore,mftypes,GCparams,MF_indices);

%% simulate the granule cells (this will take many minutes)

% we simulate in batches so as not to go nuts with memory
batch_size = 500;

% each cell is simulated this many times. there can be different outcomes
% on different runs because noise is added to the EPSPs
num_trials = 50;

% example GC so that we can make kernels
GC_model = initGCModel;

min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
tsteps      = length(min_t+GC_model.dt:GC_model.dt:max_t);

% precompute the EPSP kernels
dt                  = GC_model.dt;
win                 = 200;                  %need abt 200ms for slow*membrane
tran                = dt:dt:win;
kernel_fast         = 1/(GC_model.tau_s - GC_model.tau_m)*(exp(-tran/GC_model.tau_s) - exp(-tran/GC_model.tau_m));
kernel_slow         = 1/(GC_model.tau_s_slow - GC_model.tau_m)*(exp(-tran/GC_model.tau_s_slow) - exp(-tran/GC_model.tau_m));
spikes.kernel_fast  = kernel_fast;
spikes.kernel_slow  = kernel_slow;

% some statistics we'll keep track of
peak_vm = zeros(N,1);
spikes_per_command = zeros(N,1);
activecount = 0;

% now loop over batches and cells
for rep = 1:N/batch_size
    
    cellnames       = celltypes((rep-1)*batch_size + 1 : rep*batch_size);
    meanSpikes      = zeros(batch_size, tsteps);
    thrstore        = zeros(batch_size,1);
    vm              = cell(batch_size,1);
    sptimes         = cell(batch_size,1);
    GC_models_batch = GC_models((rep-1)*batch_size + 1 : rep*batch_size);
    
    waitbar(rep / (N/batch_size));

    % use parfor here if you want
    for ii=1:batch_size
        
        [raster, vm{ii}]          = simulate_spike_raster(GC_models_batch{ii}, rspstore, spikes, num_trials);
        sptimes{ii}               = sparse(vm{ii}>=0);
        meanSpiking               = mean(raster);
        
        if nnz(meanSpiking) ~= 0
            activecount = activecount + 1;
        end
        
        meanSpikes(ii,:)         = meanSpiking;
        thrstore(ii)             = GC_models_batch{ii}.v_thresh;
        
        peak_vm((rep-1)*batch_size + ii) = max(vm{ii}(:));
        spikes_per_command((rep-1)*batch_size + ii) = mean(meanSpiking);
    end
    
    meanSpikes = sparse(meanSpikes);
    
    % use this commented command if you want to save individual spiking
    % trials from the simulated GCs (you will need this Abby, you probably
    % won't Armen)
    save([savepth expdir '/' expname '_' num2str(rep)],'meanSpikes','thrstore','sptimes','cellnames','GCparams','GC_models_batch');
%     save([savepth expdir '/' expname '_meanonly_' num2str(rep)],'meanSpikes','thrstore','cellnames','GCparams','GC_models_batch');
end
save([savepth expdir '/pop_summary.mat'],'GC_models','activecount','celltypes','peak_vm','spikes_per_command');

