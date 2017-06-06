function [w_fast, w_slow] = get_tonic_weights_fixed(GCparams)

%pick a random ratio of fast and slow currents from our fit
%weight matrix (using all nonzero entries)
Wfast       = GCparams.Wsparse(:,:,1);
Wsum        = sum(GCparams.Wsparse,3);
ind         = find(Wsum);
Wrat        = Wfast(ind)./Wsum(ind);
Wrat        = Wrat(randi(length(Wrat),1,1));


%and pick an EPSP size:
EPSP        = 0;
while(EPSP<=0)
    EPSP    = GCparams.tonicEPSP.mean + randn(1)*GCparams.tonicEPSP.std; %distribution of EPSP peak heights (from Nate)
end

S_f         = 9.5090;
S_s         = 58.6443;

%this scaling gives an EPSP with total height EPSP, and
%Wrat % of total current injected via the slow synapse
%(I have a pdf of this derivation somewhere)
w_fast      = EPSP*S_f*S_s*Wrat/(S_s*Wrat + S_f*(1-Wrat));
w_slow      = EPSP*S_f*S_s*(1-Wrat)/(S_s*Wrat + S_f*(1-Wrat));