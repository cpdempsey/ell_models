function post_kernel = recompute_EPSP(GC_model,t,spikes)
%let's say there was a GC spike at time t; we need to recompute the synaptic
%input from that timepoint. this means take all the mf inputs within say a
%100ms window preceding the spike, and compute their remaining effect on
%the cell given what portion of their synaptic conductance presists.
%
%spikes is the original spike train, t is the time of the spike, and
%GC_model holds all the time constants.

dt          = GC_model.dt;
tau_s       = GC_model.tau_s;
tau_s_slow  = GC_model.tau_s_slow;

win         = 100; %ms. How far back to look for input history.


inds=find((spikes.sptimes<t).*(spikes.sptimes>max(t-win/dt,1))); %the spikes we care about (those recent enough to still be injecting current)


tran                = dt:dt:win;
post_kernel_fast    = zeros(size(spikes.kernel_fast));  %these are gonna store the fast and slow currents injected by previous spikes
post_kernel_slow    = zeros(size(spikes.kernel_fast));

for i=inds'
    
    %spikes.channels tells us which weights to use (ie who spiked)
    if(i<length(spikes.sptimes))
        Ws_fast = GC_model.Ws(spikes.channels(spikes.chan_id(i):spikes.chan_id(i+1)-1));
        Ws_slow = GC_model.Ws(spikes.channels(spikes.chan_id(i):spikes.chan_id(i+1)-1)+3);
    else
        Ws_fast = GC_model.Ws(spikes.channels(spikes.chan_id(i):end));
        Ws_slow = GC_model.Ws(spikes.channels(spikes.chan_id(i):end)+3);
    end
    
    %scale is how much of the EPSP is left (ie the height of the kernel
    %given a spike happened (t-spikes.spiketimes(i)) timesteps ago.)
    scale = exp(-tran(t-spikes.sptimes(i))/tau_s);
    remaining_input_fast = spikes.kernel_fast(1:win/dt)*scale;
    
    scale = exp(-tran(t-spikes.sptimes(i))/tau_s_slow);
    remaining_input_slow = spikes.kernel_slow*scale;

    %and now we add that spike's effects to the record! (it's not really
    %necessary to separate into fast and slow components here, but it makes
    %the code easier to read. Such as it is.)
    post_kernel_fast(1:length(remaining_input_fast)) = post_kernel_fast(1:length(remaining_input_fast)) ...
                        + sum(bsxfun(@times,Ws_fast',remaining_input_fast),1)*spikes.spscale(i);
                    
    post_kernel_slow(1:length(remaining_input_slow)) = post_kernel_slow(1:length(remaining_input_slow)) ...
                        + sum(bsxfun(@times,Ws_slow',remaining_input_slow),1)*spikes.spscale(i);

end

post_kernel = post_kernel_fast + post_kernel_slow; %weights were already multiplied in above