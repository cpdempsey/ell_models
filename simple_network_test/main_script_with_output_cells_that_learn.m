
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
Nmg     = 2;
Ngc     = size(gcBasis,1);
T       = size(gcBasis,2);

% weights from GCs to MG cells
Wb = zeros(Ngc,Nmg) - 0.05*randn(Ngc,Nmg);

% feedforward weights from MG cells to output cells
Wn= ones(Nmg,No);
Wn = [0.3 0.3 ; 1 1];

% feedback weights from output cells to MG cells
Wf = zeros(No,Nmg);
Wf = 1*[1  1 ; 1 1];

% feedback weights from MG cells to themselves
Wfm = zeros(Nmg,Nmg);
Wfm = 25*[1 0 ; 0 1];

% feedback weights from output cells to output cells
Wfo = zeros(No,No);
Wfo = 25*[1 0 ; 0 1];

%  weights from GCs to output cells
Wo = zeros(Ngc,No);
Wo(:,1) = 0.1*ones(Ngc,1);
Wo(:,2) = -0.1*ones(Ngc,1);

% sensory input to be cancelled
Sbase = get_amp_response(max_freq);
S = repmat(Sbase,[Nmg 1]);
S(2,:) = -S(2,:);

So = -repmat(Sbase,[No 1]);

% membrane voltage of the MG cells
M = S + Wb'*gcBasis;

% output layer activity
O = Wn'*M + So + Wo'*gcBasis;

% learning rate of MG cells
alpha = 0.000001;

% learning rate of output cells
beta = 8e-11;

% learning rate of MG cells on own voltage
gamma = 0.001;

% learning signals
L   = Wf'*O + gamma*Wfm'*M;
Lo  = Wfo'*O;


%% learning

% Wb = zeros(Ngc,Nmg);

oneT    = ones(T,Nmg);
oneTo   = ones(T,No);

clear wnorm popOutput residual residualMG residualLearning popOutput MGhist learningSig

dp = 1;
dm = 300;

Nruns = 1000;

w0 = get_rate_scales(gcBasis);
w0 = repmat(w0,[size(Wb,1) 1]);

M = S + Wb'*gcBasis;
Linit = Wf'*O;

overallLR = 4;

for ii=1:Nruns
        
    O = Wn'*M + So + Wo'*gcBasis;

    L = Wf'*O + gamma*Wfm'*M;
    
    Lo  = Wfo'*O;

    dW = dp*gcBasis*oneT - dm*gcBasis*L';
    dWo = dp*gcBasis*oneTo - dm*gcBasis*Lo';
    
    Wb = Wb + overallLR*alpha*dW;
    Wo = Wo + overallLR*beta*dWo;

    wnorm(ii) = norm(Wb);
    
    M = S + Wb'*gcBasis;
    O = Wn'*M + So + Wo'*gcBasis;
    
    popOutput(:,:,ii) = O;
    learningSig(:,ii) = L(1,:);
    MGhist(:,:,ii)  = M;
    
    residual(ii) = sum(popOutput(1,:,ii).^2) / sum(popOutput(1,:,1).^2);
    
    residualLearning(ii) = sum(learningSig(:,ii).^2) / sum(learningSig(:,1).^2);
    residualMG(ii) = mean(sum(MGhist(:,:,ii).^2,2)) / mean(sum(MGhist(:,:,1).^2,2));
    
end

L = Wf'*O;
Lo  = Wfo'*O;

negI = Wb'*gcBasis;

negIO = Wo'*gcBasis;


%%

figure;
subplot(4,2,[1 3]); hold on
ylabel('cell 1 Vm');
plot(S(1,:)); 
plot(negI(1,:));
plot(S(1,:) + negI(1,:));
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[5 7]); hold on
ylabel('cell 2 Vm');
plot(S(2,:));
plot(negI(2,:));
plot(S(2,:) + negI(2,:));
set(gca,'tickdir','out','box','off');

subplot(4,2,[2 4]); hold on
ylabel('BS rate');
plot(Linit'); 
plot(L'); 
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[6 8]);
ylabel('summed response'); hold on
Nplot = 20;
colors = winter(Nplot);
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(popOutput(1,:,plot_idx(ii)),'color',colors(ii,:));
end
plot(negIO(1,:),'color',[0 0 0],'linewidth',2);
set(gca,'tickdir','out','box','off');

setPaperPos([177 58 933 740]);

figure;
ax=plotyy(1:Nruns,residual,1:Nruns,wnorm); hold(ax(1),'on'); hold(ax(2),'on') 
plot(ax(1),1:Nruns,residualLearning,'color','green');
plot(ax(1),1:Nruns,residualMG,'color',[0.2 0.2 0.2]);
set(ax(1),'ylim',[0 1.3]);

figure; hold on
Nplot = 20;
colors = winter(Nplot);
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(learningSig(:,plot_idx(ii)),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

figure; 
subplot(1,3,1);hold on
Nplot = 20;
colors = winter(Nplot);
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(learningSig(:,plot_idx(ii)),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

subplot(1,3,2);hold on
for ii=1:Nplot
    plot(squeeze(MGhist(1,:,plot_idx(ii))),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

subplot(1,3,3);hold on
for ii=1:Nplot
    plot(squeeze(popOutput(1,:,plot_idx(ii))),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

setPaperPos([37 293 1393 357]);

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















