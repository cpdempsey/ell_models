function [new_thresh, delta_baseline] = tonic_thr_adjust(GC_model, MF_indices, spikes, rspstore)

%figure out approx firing rate of the tonic mossy, convolve with syn kernel
%to get baseline change in Vm, dance the dance of the free and the just

delta_baseline = 0;
dt=GC_model.dt;

for i=1:(length(GC_model.Ws)/2)
    if  ismember(GC_model.mf_input(i),MF_indices.tonic)
        trnum = randi(size(rspstore{GC_model.mf_input(i)},1));
        input = rspstore{GC_model.mf_input(i)}(trnum,:);
        sptrain     = find(input);
        meanISI     = mean(sptrain(2:end)-sptrain(1:end-1));  %in dt's
        if isnan(meanISI)
            continue
        end
        tonic_sc    = 1/(meanISI*dt); %this is the height of the equivalent rectangle of the spike train in sptrain
        
        fast_change = sum(spikes.kernel_fast*dt)*tonic_sc*GC_model.Ws(i);
        slow_change = sum(spikes.kernel_slow*dt)*tonic_sc*GC_model.Ws(i+(length(GC_model.Ws)/2));
        
        delta_baseline = delta_baseline + fast_change + slow_change;
          
    elseif  ismember(GC_model.mf_input(i),MF_indices.pause)
        trnum = randi(size(rspstore{GC_model.mf_input(i)},1));
        input = rspstore{GC_model.mf_input(i)}(trnum,:);
        sptrain     = find(input);
        ISIs        = diff(sptrain);
        [~,largest_interval] = max(ISIs);
        ISIs = ISIs(largest_interval+1:end);
        meanISI     = mean(ISIs);  %in dt's
        if isnan(meanISI)
            continue
        end

        pause_sc    = 1/(meanISI*dt); %this is the height of the equivalent rectangle of the spike train in sptrain
        
        fast_change = sum(spikes.kernel_fast*dt)*pause_sc*GC_model.Ws(i);
        slow_change = sum(spikes.kernel_slow*dt)*pause_sc*GC_model.Ws(i+(length(GC_model.Ws)/2));
        
        delta_baseline = delta_baseline + fast_change + slow_change;
    end
    
end

new_thresh = GC_model.v_thresh + delta_baseline;

end