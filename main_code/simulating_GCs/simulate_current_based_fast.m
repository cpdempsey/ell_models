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
spiked_early = zeros(1,length(GC_model.mf_input));
for ii=1:length(GC_model.mf_input)
    if(GC_model.mf_input(ii))
        if(rand(1)<=GC_model.MF_prob(ii))
            temp = draw_MF_input(rspstore{GC_model.mf_input(ii)});
            inputs(ii,:) = temp;
            spiked_early(ii) = nnz(temp(:,1:floor((-min_t-5)/dt)));
        end
    end
end


inputs(contains(GC_model.modeltype,'pause'),:) = ...
    circshift(inputs(contains(GC_model.modeltype,'pause'),:),[0 GC_model.jitter]);

[spikes.channels, spikes.sptimes]       = find(inputs);
[spikes.sptimes, spikes.chan_id]        = unique(spikes.sptimes);

spikes = addEPSPnoise(spikes,GC_model);       % add noise to EPSP size!

GC_spike_times=[];


% adjustments for cells that have pause or tonic input:
% if cell gets a tonic or pause input simulate some pre-command input
% to get it up to baseline correctly and avoid a period of integration
% to baseline after the command. I'm arbitrarily doing this for 0.5*tsteps
% because that's about enough time to get up to baseline for most cells
% I think, given time constants


%%%%%%%%%%%%
%%%%%%%%%%%%

% This is the new section %

if  any(spiked_early)  %any(contains(GC_model.modeltype,{'pause','tonic'}))
    GC_spike_times=[];

    GC_model.tonic = true;
    
    pre_length = 50/dt; % length of pre-command integration time to get to baseline
    
    ltemp = tsteps+pre_length;
    
    modeltrace = GC_model.El*ones(1,ltemp);
    % first find the pause and/or tonic channels
    tp_channels = find(spiked_early);%find(contains(GC_model.modeltype,{'pause','tonic'})) ;
    
    spikes_pre.sptimes = [];
    spikes_pre.channels = [];
    spikes_pre.spscale = [];
    
    for kk = 1:length(tp_channels)
        
                if contains(GC_model.modeltype{tp_channels(kk)},{'tonic'})
                    sp_int = diff(find(inputs(tp_channels(kk),:)));
                else
                    sp_int = diff(find(inputs(tp_channels(kk),1:abs(min_t/dt))));
                end
        if ~isempty(sp_int)
            
            sp_temp = sp_int(randi(length(sp_int)));
%             sp_temp = randsample(sp_int,1,true);
            
            while sp_temp(end) < pre_length
                sp_temp = [sp_temp sp_temp(end)+sp_int(randi(length(sp_int)))];
            end
            sp_temp = sp_temp(1:end-1);
        else
            sp_temp = [];
        end
        spikes_pre.sptimes  = [spikes_pre.sptimes sp_temp];
        spikes_pre.channels = [spikes_pre.channels tp_channels(kk)*ones(1,length(sp_temp))];
        
    end
    [~, Isort] = (sort(spikes_pre.sptimes));
    spikes_pre.sptimes = (spikes_pre.sptimes(Isort))';
    spikes_pre.channels = (spikes_pre.channels(Isort))';
    [spikes_pre.sptimes, spikes_pre.chan_id] = unique(spikes_pre.sptimes);    
    spikes_pre.spscale = ones(size(spikes_pre.sptimes));
    
    spikes_pre.kernel_fast = spikes.kernel_fast;
    spikes_pre.kernel_slow = spikes.kernel_slow;
    
 
i=1;
while i<=length(spikes_pre.sptimes)
    %add the next EPSP to the GC trace
    updatewin               = spikes_pre.sptimes(i)+1:min(spikes_pre.sptimes(i)+win/dt,ltemp); %bins of trace to update
    spwin                   = 1:min(ltemp-spikes_pre.sptimes(i),win/dt);                   %bins of EPSP to use for that update
    EPSP                    = compute_EPSP(GC_model,spikes_pre,i);                              %returns EPSP adjusted by apropriate weight (fast + slow components)
    modeltrace(updatewin)   = modeltrace(updatewin) + EPSP(spwin);

    %did the new EPSP cause a spike?
    if(i<length(spikes_pre.sptimes)&&(spikes_pre.sptimes(i+1)-spikes_pre.sptimes(i))<length(updatewin))
        GCspike = find(modeltrace(updatewin(1:spikes_pre.sptimes(i+1)-spikes_pre.sptimes(i)))>thr,1); %only look in window before next external spike
    else
        GCspike = find(modeltrace(updatewin)>thr,1);
    end
    
    %if it did:
    count=0;
    while((~isempty(GCspike)))
        count=count+1;        
        GCspike = GCspike + spikes_pre.sptimes(i);          %convert GCspike to an absolute time
        GC_spike_times = [GC_spike_times GCspike-1];    %add it to our record of spikes
        modeltrace(GCspike-1:end)     = GC_model.El;   %and reset the GC to El from spike time onwards

        GCspike = GCspike+GC_model.tRefrac;

        if GCspike<ltemp
            
            newinput                = recompute_EPSP(GC_model,GCspike,spikes_pre); 
            updatewin               = GCspike+1:min(GCspike+length(newinput),length(modeltrace));
            spwin                   = 1:min(length(newinput),length(modeltrace)-GCspike);
            modeltrace(updatewin)   = modeltrace(updatewin) + newinput(spwin);
               
            i = find(spikes_pre.sptimes>GCspike+1,1);   %find the next MF spike
            if(isempty(i))                          %if there are no spikes to come, things are easy
                i = length(spikes_pre.sptimes)+1;
                GCspike = find(modeltrace(spikes_pre.sptimes(end):max(updatewin))>thr,1)';
            else                                    %otherwise only look in window before next spike
                GCspike = find(modeltrace(spikes_pre.sptimes(i-1):spikes_pre.sptimes(i))>thr,1)';
            end
            
            i = i-1;
            
        else %and if we've reached the end of our simulation, empty GCspike so we can exit the while loop!
            GCspike=[];
        end
    end
    i=i+1;  
end    
    
%%


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
