function [modeltrace,spikes] = simulate_current_based_fast(GC_model,rspstore,spikes)

dt          = GC_model.dt;

min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
tsteps      = length(min_t+dt:dt:max_t);
% inputs      = zeros(length(GC_model.mf_input),tsteps+floor(tsteps/2));
win=200;  % weird number is to get a length that's a power of 2, which makes ffts faster

thr         = GC_model.v_thresh;

%pick our inputs

% select mossy fibre inputs at random from the cell's input types
inputs     = zeros(length(GC_model.mf_input),tsteps);
for ii=1:length(GC_model.mf_input)
    if(GC_model.mf_input(ii))
        if(rand(1)<=GC_model.MF_prob(ii))
            temp = draw_MF_input(rspstore{GC_model.mf_input(ii)});
            inputs(ii,:) = temp;
        end
    end
end


GC_spike_times=[];
modeltrace=GC_model.El*ones(1,tsteps*1.5);

[spikes.channels, spikes.sptimes]       = find(inputs);
[spikes.sptimes, spikes.chan_id]        = unique(spikes.sptimes);

spikes = addEPSPnoise(spikes,GC_model);       % add noise to EPSP size!

% now simulate!
i=1;
while i<=length(spikes.sptimes)
    %add the next EPSP to the GC trace
    updatewin               = spikes.sptimes(i)+1:min(spikes.sptimes(i)+win/dt,tsteps*1.5); %bins of trace to update
    spwin                   = 1:min(tsteps*1.5-spikes.sptimes(i),win/dt);                   %bins of EPSP to use for that update
    EPSP                    = compute_EPSP(GC_model,spikes,i);                              %returns EPSP adjusted by apropriate weight (fast + slow components)

    modeltrace(updatewin)   = modeltrace(updatewin) + EPSP(spwin);


    %did the new EPSP cause a spike?
    if(i<length(spikes.sptimes)&&(spikes.sptimes(i+1)-spikes.sptimes(i))<length(updatewin))
        GCspike = find(modeltrace(updatewin(1:spikes.sptimes(i+1)-spikes.sptimes(i)))>thr,1); %only look in window before next external spike
    else
        GCspike = find(modeltrace(updatewin)>thr,1);
    end
    
    %if it did:
    count=0;
    while((~isempty(GCspike)))
        count=count+1;        
        GCspike = GCspike + spikes.sptimes(i);          %convert GCspike to an absolute time
        GC_spike_times = [GC_spike_times GCspike-1];    %add it to our record of spikes
        modeltrace(GCspike-1:end)     = GC_model.El;   %and reset the GC to El from spike time onwards
        
        %we're gonna use GCspike to indicate the current time of the
        %model. nothing happens during the refractory period, so skip to
        %the end and update the GC Vm from there:
        GCspike = GCspike+GC_model.tRefrac;
        
        %if we're not done simulating, we have to compute what happens
        %after the refractory period is over:
        if(GCspike<(((GC_model.max_t-GC_model.min_t)/GC_model.dt)*1.5))
            
            %compute the effect of any recent spikes (ie those still
            %injecting current):
            newinput                = recompute_EPSP(GC_model,GCspike,spikes); 
            updatewin               = GCspike+1:min(GCspike+length(newinput),length(modeltrace));
            spwin                   = 1:min(length(newinput),length(modeltrace)-GCspike);
            modeltrace(updatewin)   = modeltrace(updatewin) + newinput(spwin);
                
                
                %see if any of those recent spikes will cause the GC to spike
            %again (this can happen if the slow current is large.) Only
            %check in the time window before the next MF spike, as
            %that will change subsequent spiketimes:
            i = find(spikes.sptimes>GCspike+1,1);   %find the next MF spike
            if(isempty(i))                          %if there are no spikes to come, things are easy
                i = length(spikes.sptimes)+1;
                GCspike = find(modeltrace(spikes.sptimes(end):max(updatewin))>thr,1)';
            else                                    %otherwise only look in window before next spike
                GCspike = find(modeltrace(spikes.sptimes(i-1):spikes.sptimes(i))>thr,1)';
            end
            
            %and with that we've marked the time of the next evoked spike
            %(by storing it in GCspike.) The next iteration of the while
            %loop will implement its effects and carry on. Since i indexes
            %the *next* MF spike (after GCspike) and we're still
            %working out the effects of the most recent spike, subtract 1
            %from i (so the next iteration of the while loop is still
            %working from the most recent spike.)
            i = i-1;
            
        else %and if we've reached the end of our simulation, empty GCspike so we can exit the while loop!
            GCspike=[];
        end
    end
    i=i+1;  %hooray, we've found all the GC spikes happening between MF spikes i and i+1.
            %now let's move on to the next interval.
end

modeltrace(GC_spike_times)=0; %make the spikes look like spikes
modeltrace=modeltrace(1:tsteps); % and cut off the half-interval we added to the front of our simulation

end
