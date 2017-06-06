function post_kernel = compute_EPSP(GC_model,spikes,i)

%channels tells us which weights to use (ie who spiked)
if(i<length(spikes.sptimes))
    Ws_fast = GC_model.Ws(spikes.channels(spikes.chan_id(i):spikes.chan_id(i+1)-1));
    Ws_slow = GC_model.Ws(spikes.channels(spikes.chan_id(i):spikes.chan_id(i+1)-1)+3);
else
    Ws_fast = GC_model.Ws(spikes.channels(spikes.chan_id(i):end));
    Ws_slow = GC_model.Ws(spikes.channels(spikes.chan_id(i):end)+3);
end


scaled_kernel_fast = sum(bsxfun(@times,Ws_fast',spikes.kernel_fast),1)*spikes.spscale(i); %spscale is the EPSP noise
scaled_kernel_slow = sum(bsxfun(@times,Ws_slow',spikes.kernel_slow),1)*spikes.spscale(i);

post_kernel = scaled_kernel_fast + scaled_kernel_slow; %weights were already multiplied in above

end