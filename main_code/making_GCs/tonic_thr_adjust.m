function [new_thresh, delta_baseline] = tonic_thr_adjust(GC_model, MF_indices, spikes, rspstore)

%figure out approx firing rate of the tonic mossy, convolve with syn kernel
%to get baseline change in Vm, dance the dance of the free and the just

delta_baseline = 0;
dt=GC_model.dt;

for i=1:(length(GC_model.Ws)/2)
    if(sum(MF_indices.tonic==GC_model.mf_input(i)))
        sptrain     = find(draw_MF_input(rspstore{GC_model.mf_input(i)}));
        meanISI     = mean(sptrain(2:end)-sptrain(1:end-1));  %in dt's
        
        tonic_sc    = 1/(meanISI*dt); %this is the height of the equivalent rectangle of the spike train in sptrain
        
        fast_change = sum(spikes.kernel_fast*dt)*tonic_sc*GC_model.Ws(i);
        slow_change = sum(spikes.kernel_slow*dt)*tonic_sc*GC_model.Ws(i+(length(GC_model.Ws)/2));
        
        delta_baseline = delta_baseline + fast_change + slow_change;
    end
end

new_thresh = GC_model.v_thresh + delta_baseline;

end