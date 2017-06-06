function initWeightsRealistic(gcBasis,Rs,Rb,Wn,Wf,El,dm)

Nmg = size(Wn,1);
No = size(Wn,2);
T = size(gcBasis,2);

% weights from GCs to MG cells
Wb = zeros(Ngc,Nmg) + 0.02*rand(Ngc,Nmg);

% membrane voltage of the MG cells
Vinit = El + Wb'*gcBasis;

% output rate
Vo = El + Wn'*Rs(Vinit);
Rs_out = Rs(Vo);

oneT    = ones(T,Nmg);

Nruns = 5000;

kr = 0.001;

fb = Wf'*Rs_out;

for ii=1:Nruns
    
    waitbar(ii/Nruns);
    
    V = El + Wb'*gcBasis;
    
    Rs_current =    Rs(V);
    Rb_current =    Rb(V,fb);
    
    Vo = El + Wn'*Rs_current;
    Rs_out = Rs(Vo);
    
    fb = Wf'*Rs_out;
    
    dW = gcBasis*oneT - dm*gcBasis*Rb_current';
    %     dWo = dp*gcBasis*oneTo - dm*gcBasis*Lo';
    
    Wb = Wb + kr*dW;
    %     Wo = Wo + overallLR*beta*dWo;
    
    % bounce off boundary
    Wb(Wb<0) = -Wb(Wb<0);
    

end