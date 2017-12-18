function [GC_models, celltypes] = generateGCPopAnn(N,rspstore,mftypes,GCparams,MF_indices)

% generate distance to threshold
thr_mean = 17.5;
% thr_mean = 13.5;
thr_var = 6.8^2 + 4.9^2;

GC_models       = cell(N,1);

celltypes = cell(N,1);

% make a standard GC model
GC_model = initGCModel;

%precompute the EPSP kernels!
dt                  = GC_model.dt;
win                 = 200;                  %need abt 200ms for slow*membrane
tran                = dt:dt:win;
kernel_fast         = 1/(GC_model.tau_s - GC_model.tau_m)*(exp(-tran/GC_model.tau_s) - exp(-tran/GC_model.tau_m));
kernel_slow         = 1/(GC_model.tau_s_slow - GC_model.tau_m)*(exp(-tran/GC_model.tau_s_slow) - exp(-tran/GC_model.tau_m));
spikes.kernel_fast  = kernel_fast;
spikes.kernel_slow  = kernel_slow;

for ii = 1:N
    
    GC_models{ii} = initGCModel;
      
    %%%% here we draw input types and weights using parameters fit by Ann
    [MF_input, Ws, cellname, GC_models{ii}.modeltype] = getInputMix(GCparams,MF_indices);
        
    celltypes{ii} = cellname;

    GC_models{ii}.cellname = cellname;
    GC_models{ii}.mf_input = MF_input;
    
    GC_models{ii}.Ws = Ws; 
                  
    GC_models{ii}.v_thresh = GC_models{ii}.El;
    while(GC_models{ii}.v_thresh<=GC_models{ii}.El+3) %don't let threshold go below El + 3
        GC_models{ii}.v_thresh = GC_models{ii}.El + thr_mean + randn(1)*sqrt(thr_var);
    end
    
    tonic_on = any(cellfun(@(x) ~isempty(x) , strfind(mftypes(GC_models{ii}.mf_input(GC_models{ii}.mf_input~=0)),'tonic')));
    pause_on = any(cellfun(@(x) ~isempty(x) , strfind(mftypes(GC_models{ii}.mf_input(GC_models{ii}.mf_input~=0)),'pause')));
    
    if tonic_on || pause_on
        [GC_models{ii}.v_thresh,d_baseline] = tonic_thr_adjust(GC_models{ii}, MF_indices, spikes, rspstore);
        GC_models{ii}.baseline = GC_models{ii}.El + d_baseline;
    else
        GC_models{ii}.baseline = GC_models{ii}.El;
    end

end

end