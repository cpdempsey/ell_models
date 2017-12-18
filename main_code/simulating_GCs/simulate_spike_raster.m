function [raster, vm] = simulate_spike_raster(GC_model,rspstore, spikes, numtrials)

dt          = GC_model.dt;

min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
tsteps      = length(min_t+dt:dt:max_t);


raster=zeros(numtrials,tsteps);
vm=zeros(numtrials,tsteps);

jitter_mean = 5; %ms
jitter_std = 5;

GC_model.jitter = round( (1/dt)*(jitter_mean + jitter_std*randn(1,1)) );

for i=1:numtrials
    [modeltrace] = simulate_current_based_fast(GC_model,rspstore, spikes);
    
%     if any(contains(GC_model.modeltype,'pause'))
%         disp('pause');
%     end
    
    sptrain = zeros(1,length(modeltrace));
    sptrain(modeltrace==0) = 1;
    raster(i,:) = sptrain;
    vm(i,:)     = modeltrace;
end

end