function spikes = addEPSPnoise(spikes,GC_model)

%our spikes get convolved with the membrane (spike->spike/GC_model.Wscale),
%and then multiplied by the synaptic weight-- so we include these two terms
%in the variance addition so that our final EPSP has variance
%GC_model.EPSPvar.

spikes.spscale = zeros(size(spikes.sptimes));


varscale=zeros(1,length(GC_model.Ws));
for Wnum=1:(length(GC_model.Ws)/2)
    if(GC_model.Ws(Wnum))
        varscale(Wnum) = varscale(Wnum) + GC_model.Ws(Wnum)/GC_model.Wscale_fast;
    end
    if(GC_model.Ws(Wnum+3))
        varscale(Wnum) = varscale(Wnum) + GC_model.Ws(Wnum+3)/GC_model.Wscale_slow;
    end
    if varscale(Wnum) == 0
        varscale(Wnum) = 1;
    end
end

for i=1:length(spikes.sptimes)
    noise=-inf;
    if(i<length(spikes.sptimes))
        for j=spikes.chan_id(i):spikes.chan_id(i+1)-1
            while((noise+1)<0) %don't go negative!
                noise = randn(1)*sqrt(GC_model.EPSPvar/ varscale(spikes.channels(j))^2);
            end
        end
    else
        for j=spikes.chan_id(i):length(spikes.channels)
            while((noise+1)<0) %don't go negative!
                noise = randn(1)*sqrt(GC_model.EPSPvar/ varscale(spikes.channels(j))^2);
            end
        end
    end
    spikes.spscale(i) = 1 + noise;
end

end


