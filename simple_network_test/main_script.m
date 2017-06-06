
% simulate MG cell learning

% plotting stuff:
set(0,'DefaultAxesLinewidth',1.5)
set(0,'DefaultLineLinewidth',3)
set(0,'defaultaxesfontsize',24);

rmpath(genpath('~/Documents/MATLAB/'));
cd('/Users/seigen/Documents/MATLAB/ell_models');
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
F = [0.2 0.9; 
    -0.2 0.2];

% sensory input to be cancelled
S = get_amp_response(max_freq);
S = repmat(S,[Nmg 1]);

% membrane voltage of the MG cells
Vmg = S' + gcBasis'*W;


%% learning

W = zeros(Ngc,Nmg);
Vmg = S' + gcBasis'*W;

oneT = ones(T,Nmg);

clear wnorm popOutput residual

dp = 1;
dm = 300;
alpha = 0.0000001;

Nruns = 5000;

w0 = get_rate_scales(gcBasis);
w0 = repmat(w0,[size(W,1) 1]);

Vmg = S' + gcBasis'*W;
Linit = Vmg*F;

for ii=1:Nruns
        
    L = Vmg*F;

    dW = dp*gcBasis*oneT - dm*(gcBasis*L);
    
    W = W + alpha*dW;
    
    wnorm(ii) = norm(W);
    
    Vmg = S' + gcBasis'*W;
    
    popOutput(:,ii) = mean(Vmg,2);
    
    residual(ii) = sum(popOutput(:,ii).^2) / sum(popOutput(:,1).^2);
    
end

L = Vmg*F;
negI = gcBasis'*W;

%%
figure;
subplot(4,2,[1 3]); hold on
ylabel('cell 1 Vm');
plot(S(1,:)); 
plot(negI(:,1));
plot(S(1,:)' + negI(:,1));
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[5 7]); hold on
ylabel('cell 2 Vm');
plot(S(2,:));
plot(negI(:,2));
plot(S(2,:)' + negI(:,2));
set(gca,'tickdir','out','box','off');

subplot(4,2,[2 4]); hold on
ylabel('BS rate');
plot(Linit); 
plot(L); 
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[6 8]);
ylabel('summed response'); hold on
Nplot = 20;
colors = winter(Nplot);
plot(mean(S)); 
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(popOutput(:,plot_idx(ii)),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

setPaperPos([177 58 933 740]);

figure;
plotyy(1:Nruns,residual,1:Nruns,wnorm); 

%% look at response to some new signal post learning

Snew = -S;
Vmg = Snew' + gcBasis'*W;

L = Vmg*F;
negI = gcBasis'*W;

figure;
subplot(4,2,[1 3])
plot(Snew(1,:)); hold on
plot(negI(:,1));
plot(Snew(1,:)' + negI(:,1));
set(gca,'tickdir','out','box','off');

subplot(4,2,[5 7]);
plot(Snew(2,:)); hold on
plot(negI(:,2));
plot(Snew(2,:)' + negI(:,2));
set(gca,'tickdir','out','box','off');

subplot(4,2,[2 4]);
plot(L); 
set(gca,'tickdir','out','box','off');

subplot(4,2,[6 8]);
plot(mean(Snew)); hold on
plot(mean(Vmg,2))

setPaperPos([177 58 933 740]);















