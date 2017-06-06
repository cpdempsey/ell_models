
% simulate MG cell learning

% plotting stuff:
set(0,'DefaultAxesLinewidth',1.5)
set(0,'DefaultLineLinewidth',3)
set(0,'defaultaxesfontsize',24);

rmpath(genpath('~/Documents/MATLAB/'));
cd('/Users/seigen/Documents/MATLAB/ELL_network_project');
addpath(genpath('./'));

%% get a summary of the GC population activity

% parameters of kernel that we convole the GC spikes with to turn them into
% a voltage
syn_kernel      = struct('tau_on', 0.005, 'tau_off', 0.01,'dt',5e-05);

% where are the granule cells we want to use
pth = './generated_GCs/';

% load the basis
[~, gcBasis, ~, ~, summedSpikes, celltypes_firing, max_freq] = getPopulationResponse(pth,syn_kernel);

gcBasis = gcBasis';

%% delay line basis

c_w_syn    = @(x) convolve_with_synaptic_kernel(x', syn_kernel.dt,syn_kernel.tau_on, syn_kernel.tau_off);

delay_line = eye(4500);

delay_line = c_w_syn(delay_line);

delay_line = freqLimit(delay_line,0.05/1000);

gc_rates   = delay_line;


%% basic parameters

rectNonLin = @(x) max(x,0);

No      = 2;
Nmg     = 4;
Ngc     = size(gcBasis,1);
T       = size(gcBasis,2);

% E or I type for each MG cell
MGtypes = 2*(round(rand(Nmg,1))-0.5);

% weights from GCs to MG cells (excitatory)
Wb = zeros(Ngc,Nmg) + 0.05*rand(Ngc,Nmg);
w0 = 0.1*get_rate_scales(gcBasis')';
Wb = repmat(w0,[1 size(Wb,2)]);

% feedforward weights from MG cells to output cells (inhibitory)
Wn= 0.01*randn(Nmg,No);

% feedback weights from output cells to MG cells
Wf = randn(No,Nmg);

% sensory input to be cancelled
baselineSensory = 0;
scaleSensory = 0;

Sbase = get_amp_response(max_freq);
S = repmat(Sbase,[Nmg 1]);
S = bsxfun(@times,MGtypes,S);
S = scaleSensory*(S + baselineSensory);

So = repmat(Sbase,[No 1]);
So = scaleSensory*(So + baselineSensory);

% membrane voltage of the MG cells
M = S + Wb'*gcBasis;

% output layer activity
O = Wn'*M + So;

% learning signals
L   = Wf'*O;


%% learning

% Wb = zeros(Ngc,Nmg);

oneT    = ones(T,Nmg);

clear wnorm popOutput residual residualMG residualLearning Mhist learningSig

dm = 0.5;
kr = 0.0000003;

Nruns = 1000;

w0 = get_rate_scales(gcBasis);
w0 = repmat(w0,[size(Wb,1) 1]);

M = S + Wb'*gcBasis;
Linit = Wf'*O;

learningSig = zeros(size(Linit,1),size(Linit,2),Nruns);
Mhist       = zeros(size(M,1),size(M,2),Nruns);
popOutput   = zeros(size(O,1),size(O,2),Nruns);

tic
for ii=1:Nruns
    
    waitbar(ii/Nruns);
        
    O = Wn'*M + So;

    L = Wf'*O;
    
    dW = kr*(gcBasis*oneT - dm*gcBasis*L');
    
    Wb = Wb + dW;
    
%     Wb(Wb < 0) = -Wb(Wb < 0);

    wnorm(ii) = norm(Wb);
    
    M = S + Wb'*gcBasis;
    O = Wn'*M + So;
    
    popOutput(:,:,ii) = O;
    learningSig(:,:,ii) = L;
    Mhist(:,:,ii)  = M;
    
    residual(ii) = mean(var(popOutput(:,:,ii),[],2) ./ var(popOutput(:,:,1),[],2) );
    
    residualLearning(ii) = mean(var(learningSig(:,:,ii),[],2) ./ var(learningSig(:,:,1),[],2) );
    residualMG(ii) = mean(var(Mhist(:,:,ii),[],2) ./ var(Mhist(:,:,1),[],2));
    
end
time_taken = toc;

L = Wf'*O;

negI = Wb'*gcBasis;


%%


plot(residual); hold on
plot(residualLearning);








