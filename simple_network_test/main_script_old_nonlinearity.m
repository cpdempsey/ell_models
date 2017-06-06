
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

Nmg  = 2;
No   = 1;
Ngc = size(gcBasis,1);
T = size(gcBasis,2);

% weights from GCs to MG cells
W = zeros(Ngc,Nmg);

% feedback weights
% F = 0.01*randn(Nmg, Nmg);
% F = eye(Nmg) + 0.01*randn(Nmg,Nmg);
F = 1/Nmg * ones(Nmg,Nmg);

% sensory input to be cancelled
S = get_amp_response(max_freq);
S = repmat(S,[Nmg 1]);

% MG cell parameters
gamma = 0.1; % units of 1/mV
theta = -60*ones(Nmg,1);  % units of mV
theta = repmat(theta,[1, T]);
Rm = 1/30;      % max rate in units of Hz
El = -60;

% membrane voltage of the MG cells
Vmg = S + W*gcBasis' + El;

%%

RofV = @(V) 2*Rm ./ (1 + exp(-gamma*(F*V-theta)));

tt=linspace(-25,200,T);

Rb = RofV(Vmg);

plotyy(tt,Rb(1,:),tt,Vmg(1,:));

%% learning

W = zeros(Nmg,Ngc);
Vmg = S + W*gcBasis' + El;

clear wnorm

dp = 1;
dm = 30;
alpha = 0.01;

gcSum = repmat(sum(gcBasis),[Nmg 1]);

Nruns = 500;

w0 = get_rate_scales(gcBasis);
w0 = repmat(w0,[size(W,1) 1]);

for ii=1:Nruns
    
    Vmg = S + W*gcBasis' + El;

    dW = dp*gcSum - dm*(RofV(Vmg)*gcBasis);
    
    W = W + alpha*dW;
    
    wnorm(ii) = norm(Vmg);
    
end


figure;
plot(S(2,:)+El); hold on
plot(W(2,:)*gcBasis'+El);
plot(S(2,:) + W(2,:)*gcBasis'+El);





















