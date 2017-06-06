
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

No      = 1;
Nmg     = 1;
Ngc     = size(gcBasis,1);
T       = size(gcBasis,2);

% E or I type for each MG cell
MGtypes = 2*(round(rand(Nmg,1))-0.5);

% weights from GCs to MG cells (excitatory)
Wb = zeros(Ngc,Nmg) + 0.05*rand(Ngc,Nmg);
w0 = 0.1*get_rate_scales(gcBasis')';
Wb = repmat(w0,[1 size(Wb,2)]);

% feedforward weights from MG cells to output cells (inhibitory)
Wn= -0.01*rand(Nmg,No);
% Wn = [0.3 0.3 ; 1 1];

% feedback weights from output cells to MG cells
% Wf = rand(No,Nmg);
Wf = 1;
% Wf = 1*[1  1 ; 1 1];

% feedback weights from MG cells to themselves
Wfm = eye(Nmg,Nmg);
% Wfm = 25*[1 0 ; 0 1];

% feedback weights from output cells to output cells
Wfo = eye(No,No);
% Wfo = 25*[1 0 ; 0 1];

%  weights from GCs to output cells
Wo = 0.01*rand(Ngc,No);
% Wo(:,1) = 0.1*ones(Ngc,1);
% Wo(:,2) = -0.1*ones(Ngc,1);

% sensory input to be cancelled
baselineSensory = 10;

Sbase = get_amp_response(max_freq);
S = repmat(Sbase,[Nmg 1]);
S = bsxfun(@times,MGtypes,S);
S = S + baselineSensory;

So = repmat(Sbase,[No 1]);
So = So + baselineSensory;

% membrane voltage of the MG cells
M = S + Wb'*gcBasis;

% output layer activity
O = Wn'*M + So + Wo'*gcBasis;

% learning rate of MG cells
alpha = 0.000001;

% learning rate of output cells
beta = 0;

% learning rate of MG cells on own voltage
gamma = 0;

% learning signals
L   = Wf'*O + gamma*Wfm'*M;
Lo  = Wfo'*O;


%% learning

% Wb = zeros(Ngc,Nmg);

oneT    = ones(T,Nmg);
oneTo   = ones(T,No);

clear wnorm popOutput residual residualMG residualLearning Mhist learningSig

dp = 1;
dm = 0.12;

Nruns = 1000;

w0 = get_rate_scales(gcBasis);
w0 = repmat(w0,[size(Wb,1) 1]);

M = S + Wb'*gcBasis;
Linit = Wf'*O;

learningSig = zeros(size(Linit,1),size(Linit,2),Nruns);
Mhist       = zeros(size(M,1),size(M,2),Nruns);
popOutput   = zeros(size(O,1),size(O,2),Nruns);

overallLR = 100;
tic
for ii=1:Nruns
        
    O = Wn'*M + So + Wo'*gcBasis;

    L = Wf'*O + gamma*Wfm'*M;
    
    Lo  = Wfo'*O;

    dW = dp*gcBasis*oneT - dm*gcBasis*L';
    dWo = dp*gcBasis*oneTo - dm*gcBasis*Lo';
    
    Wb = Wb + overallLR*alpha*dW;
    Wo = Wo + overallLR*beta*dWo;
    
    Wb(Wb < 0) = -Wb(Wb < 0);
    Wo(Wo < 0) = -Wo(Wo < 0);

    wnorm(ii) = norm(Wb);
    
    M = S + Wb'*gcBasis;
    O = Wn'*M + So + Wo'*gcBasis;
    
    popOutput(:,:,ii) = O;
    learningSig(:,:,ii) = L;
    Mhist(:,:,ii)  = M;
    
    residual(ii) = mean(var(popOutput(:,:,ii),[],2) ./ var(popOutput(:,:,1),[],2) );
    
    residualLearning(ii) = mean(var(learningSig(:,:,ii),[],2) ./ var(learningSig(:,:,1),[],2) );
    residualMG(ii) = mean(var(Mhist(:,:,ii),[],2) ./ var(Mhist(:,:,1),[],2));
    
end
time_taken = toc;

L = Wf'*O;
Lo  = Wfo'*O;

negI = Wb'*gcBasis;

negIO = Wo'*gcBasis;


%%

figure;
subplot(4,2,[1 3]); hold on
title('MG cell 1')
plot(S(1,:)); 
plot(negI(1,:));
plot(S(1,:) + negI(1,:));
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[5 7]); hold on
title('MG cell 2')
plot(S(min(size(S,1),2),:));
plot(negI(min(size(S,1),2),:));
plot(S(min(size(S,1),2),:) + negI(min(size(S,1),2),:));
set(gca,'tickdir','out','box','off');

subplot(4,2,[2 4]); hold on
title('BS rate');
plot(Linit'); 
plot(L'); 
set(gca,'tickdir','out','box','off','xtick',[]);

subplot(4,2,[6 8]);
title('example output'); hold on
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
title('residuals');
plot(ax(1),1:Nruns,residualLearning,'color','green');
plot(ax(1),1:Nruns,residualMG,'color',[0.2 0.2 0.2]);
set(ax(1),'ylim',[0 1.3]);

figure; hold on
title('example learning signals');
Nplot = 20;
colors = winter(Nplot);
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(learningSig(1,:,plot_idx(ii)),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

figure; 
subplot(1,3,1);hold on
title('learning signals');
Nplot = 20;
colors = winter(Nplot);
plot_idx = round(linspace(1,Nruns,Nplot));
for ii=1:Nplot
    plot(learningSig(1,:,plot_idx(ii)),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

subplot(1,3,2);hold on
title('MG response e.g.');
for ii=1:Nplot
    plot(squeeze(Mhist(1,:,plot_idx(ii))),'color',colors(ii,:));
end
set(gca,'tickdir','out','box','off');

subplot(1,3,3);hold on
title('output response e.g.');
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















