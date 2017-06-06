
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

No      = 2;
Nmg     = 50;
Ngc     = size(gcBasis,1);
T       = size(gcBasis,2);

% baseline voltage
El = -67;

% max broad/simple spike rates
Rb_max = 10; % Hz
Rs_max = 50;

% steepness of voltage to rate curves
alpha_s = 0.4;
alpha_b = 1;
alpha_f = 3;

% thresholds for broad and simple spikes
theta_s = -62;
theta_b = -62;

% transfer functions
Rs =  @(V)  Rs_max ./ (1 + exp(-alpha_s*(V - theta_s)));
Rb =  @(V,fb)  Rb_max ./ (1 + exp(-alpha_b*(V - theta_b)).*exp(-alpha_f*fb));

% E or I type for each MG cell
MGtypes = 2*(round(rand(Nmg,1))-0.5);
% MGtypes = ones(size(MGtypes));

% weights from GCs to MG cells
Wb = zeros(Ngc,Nmg) + 0.02*rand(Ngc,Nmg);

% feedforward weights from MG cells to output cells
Wn = -0.01*(1/sqrt(Nmg))*rand(Nmg,No);
% Wn = [0.3 0.3 ; 1 1];

% feedback weights from output cells to MG cells
Wf = 3*(1/sqrt(No))*(1/Rs_max)*randn(No,Nmg) + 0/(sqrt(No)*Rs_max);

% sensory input to be cancelled
max_sensory = 1; %in mV
Sbase = get_amp_response(max_freq);
Sbase = max_sensory*Sbase/max(Sbase);
S = repmat(Sbase,[Nmg 1]);
S = bsxfun(@times,MGtypes,S);

So = repmat(Sbase,[No 1]);

% membrane voltage of the MG cells
Vinit = El + S + Wb'*gcBasis;

% output layer activity
O = Wn'*Vinit + So + Wo'*gcBasis;

% learning signals
L   = Wf'*O + gamma*Wfm'*Vinit;

% output rate
Vo = El + So + Wn'*Rs(Vinit);
Rs_out = Rs(Vo);

% plot transfer functions
Vplot = linspace(-70,-60,1000);
% Rbtarg = Rb_target(Rs_out);
plot(Vplot,Rb(Vplot,1)); hold on
plot(Vplot,Rs(Vplot));

dm = 0.5;

Rb_eq = 1/dm;

Veq = theta_b + log((Rb_max - Rb_eq - 1)^(-1/alpha_b));

Rs_eq_special = Rs_max / (Rb_max - Rb_eq);
Rs_eq = Rs(Veq);

%% initialise weights

Wbinit = initWeightsRealistic(gcBasis,Rs,Rb,Wn,Wf,El,dm,0.004,1500);

%% command only responses pre

Vpre = El + Wbinit'*gcBasis;

Rs_pre =    Rs(Vpre);

Vo_pre = El + Wn'*Rs_pre;
Rs_out_pre = Rs(Vo_pre);

fb_pre = Wf'*Rs_out_pre;

Rb_pre =    Rb(Vpre,fb_pre);

figure;
subplot(131); plot(Rs_out_pre');
subplot(132); plot(Rs_pre');
subplot(133); plot(Rb_pre');

setPaperPos([47 467 1286 339]);

%% sensory responses pre

V = El + Wbinit'*gcBasis + S;

Rsp =  Rs(V);

Vop = El + Wn'*Rsp + So;
Rsop = Rs(Vop);

fbp = Wf'*Rsop;

Rbp =    Rb(V,fbp);

figure;
subplot(131); plot(Rsop');
subplot(132); plot(Rsp');
subplot(133); plot(Rbp');

setPaperPos([44 29 1286 339]);

%% learning

% Wb = zeros(Ngc,Nmg);

oneT    = ones(T,Nmg);
oneTo   = ones(T,No);

clear wnorm Rshist Rbhist Rohist

Nruns = 5000;

kr = 0.003;

fb = Wf'*Rs_out;

Wb = Wbinit;

tic
for ii=1:Nruns
    
    waitbar(ii/Nruns);
    
    V = El + S + Wb'*gcBasis;
    
    Rs_current =    Rs(V);
    Rb_current =    Rb(V,fb);
    
    Vo = El + So + Wn'*Rs_current;
    Rs_out = Rs(Vo);
    
    fb = Wf'*Rs_out;
    
    dW = gcBasis*oneT - dm*gcBasis*Rb_current';
    
    Wb = Wb + kr*dW;
    
    % bounce off boundary
    Wb(Wb<0) = -Wb(Wb<0);
    
    wnorm(ii) = norm(Wb);

    Rbhist(:,:,ii)  = Rb_current;
    Rshist(:,:,ii)  = Rs_current;
    
    % compute command only responses
    Vn = El + Wb'*gcBasis;
    Rsn =    Rs(Vn);
    Von = El + Wn'*Rsn;
    Ron = Rs(Von);
    fbn = Wf'*Ron;
    Rbn =    Rb(Vn,fbn);
    
    negBhist(:,:,ii) = Rbn;
    negShist(:,:,ii) = Rsn;

    Vmghist(:,:,ii) = V;
    
    Rohist(:,:,ii)  = Rs_out;
    
    residualSimple(ii) = mean(var(Rshist(:,:,ii),[],2) ./ var(Rshist(:,:,1),[],2) );
    
    residualBroad(ii) = mean(var(Rbhist(:,:,ii),[],2) ./ var(Rbhist(:,:,1),[],2) );
    residualOut(ii) = mean(var(Rohist(:,:,ii),[],2) ./ var(Rohist(:,:,1),[],2));
    
end
time_taken = toc;

% L = Wf'*O;
% Lo  = Wfo'*O;

% negI = Wb'*gcBasis;

% negIO = Wo'*gcBasis;


%% plot a variety of MG cells start and end

figure;
count=1;
idx_plot = randsample(Nmg,15);

for ii=idx_plot'
    subplot(3,5,count)
    plot(Rshist(ii,:,1)'); hold on
    plot(Rshist(ii,:,end)');
    count = count + 1;
end

setPaperPos([45 1 1391 805]);


%% plot a variety of output cells start and end

figure;
for ii=1:No
    subplot(1,No,ii)
    plot(Rohist(ii,:,1)'); hold on
    plot(Rohist(ii,:,end)');
end

setPaperPos([39 294 1306 384]);

%% plot residuals

plot(residualSimple); hold on
plot(residualBroad)
plot(residualOut);
hl=legend({'simple','broad','output'});
set(hl,'box','off');

%% plot command only responses early and post command only


figure;
count=1;
idx_plot = randsample(Nmg,15);

for ii=idx_plot'
    subplot(3,5,count)
    plot(negBhist(ii,:,1)'); hold on
    plot(negBhist(ii,:,end)');
    count = count + 1;
end

setPaperPos([45 1 1391 805]);






