function [raster, vm] = simulate_spike_raster(GC_model,rspstore, spikes, numtrials)

dt          = GC_model.dt;

min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
tsteps      = length(min_t+dt:dt:max_t);


raster=zeros(numtrials,tsteps);
vm=zeros(numtrials,tsteps);

for i=1:numtrials
    [modeltrace] = simulate_current_based_fast(GC_model,rspstore, spikes);
    raster(i,modeltrace==0) = 1;
    vm(i,:)     = modeltrace;
end

end