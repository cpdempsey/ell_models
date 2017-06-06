% simulate MG cell learning

% plotting stuff:
set(0,'DefaultAxesLinewidth',1.5)
set(0,'DefaultLineLinewidth',3)
set(0,'defaultaxesfontsize',24);

%%

rmpath(genpath('~/Documents/MATLAB'));
cd('~/Documents/MATLAB/ELL_network_project/');
addpath(genpath('./'));


%% get a summary of the GC population activity

% parameters of kernel that we convole the GC spikes with to turn them into
% a voltage
syn_kernel      = struct('tau_on', 0.005, 'tau_off', 0.01,'dt',5e-05);

% where are the granule cells we want to use
pth = './generated_GCs/';

% load the basis
[~, gc_rates, ~, ~, summedSpikes, ~, max_freq] = getPopulationResponse(pth,syn_kernel);

%% do learnnig

% max MG response to the ampullary input in mV
peak_vm = 8;

% load ampullary response 
amp_response = get_amp_response(max_freq);
amp_response = peak_vm * amp_response / max(amp_response);

tvec = linspace(-0.025,0.2,length(amp_response));

% regularization parameter (strength of the regularization)
lambda = 0.01;

% iterations of the learning
Nruns = 5000;

% ratio between associative and non-associative plasticity. Armen, just
% turn deltaM to 0 to simulate the loss of associative plasticity. The
% result of that experiment will depend on the bounds on the weights,
% whether or not you use the scaling of the learning rate described below
% (that was used in Ann's paper), and obviously what you think the rate of
% learning is and how long the experiment goes on for.
deltaP = 1;
deltaM = -70;

% scale of the learning rule
alpha = 0.001;

dm = deltaM / size(gc_rates,1);
dp = deltaP / size(gc_rates,1);

% in Ann's paper she scales the learning rate of a cell roughly in inverse
% proportion to the amount of GC activity at the time of its peak firing,
% you can try it with and without that
wscale = get_rate_scales(gc_rates);
scale_on =  false;

if scale_on
    S = diag(wscale);
else
    S = eye(size(wscale,1));
end

% set up the various matrices used for simulating the learning
% feel free to ask me (Conor) about this if you want it explained
R = gc_rates';
N = size(R,1);
T = size(R,2);

A = vec(amp_response);
Q = lambda*eye(N) - dm*S*(R*R');
id_T = ones(T,1);

% target weight vector towards which the regularization will push you
wtarget = 0*ones(size(gc_rates,2),1);

% initial weights (this doesn't matter unless there's no regularization at
% all, i.e. unless lambda is 0)
winit   = 0*ones(size(gc_rates,2),1);

b = dp*S*R*id_T + lambda*wtarget + dm*S*R*A;

% we can calculate what the final weight vector would be after infinite
% learning with unbounded weights. If there's no regularization or bounding 
% of the weights then learning will create a perfect negative image. This
% ignores the effects of trial over trial noise in the sensory response and
% the granule cell responses. These noise levels are also an important part
% of what will really determine when learning stops, but are not explicitly
% considered here. But this might be important when thinking about learning
% rates and the level of cancellation eventually achieved.
w_equil = Q\b;
neg_i_final = gc_rates*w_equil;
canc_final = A + neg_i_final;

% GC basis elements scaled by their equilibrium weight
contribs_equil = bsxfun(@times,gc_rates,w_equil');

figure;
plot(tvec,A); hold on ; plot(tvec,canc_final); plot(tvec,neg_i_final); 
set(gca,'xlim',[tvec(1) tvec(end)],'box','off','tickdir','out');
xlabel('time (s)'); ylabel('MG Vm');

% upper and lower weight bounds in mV (this bounds things so that the lare
bd = 0.4;
ub = bd*ones(size(winit)); lb = -bd*ones(size(winit));
ub = ub / max(gc_rates(:)); lb = lb / max(gc_rates(:));

% initial weights
w = winit;

% do the learning
for ii = 1:Nruns
    
    waitbar(ii/Nruns);
    
    % update weights
    w = w - alpha*(Q*w - b);
    % we can either bounce off the bounds or stick to them (uncomment one)
    
    % bounce:
%     w(w>ub) = w(w>ub) - 2*(w(w>ub) - ub(w>ub));
%     w(w<lb) = w(w<lb) + 2*(lb(w<lb) - w(w<lb));
     
    % stick: 
    w(w>ub) = ub(w>ub);
    w(w<lb) = lb(w<lb);
    
    whist(ii,:) = w;
    
    % current negative image
    neg_i_store(ii,:) = gc_rates*w;
    canc_store(ii,:) = A + neg_i_store(ii,:)';
    
    mse_ratio(ii) = sum(canc_store(ii,:).^2) / sum(A.^2);
end

% GC basis elements scaled by their final weight
contribs = bsxfun(@times,gc_rates,w');

% plot progression of learning
idx_plot = unique(round(logspace(log(1),log(Nruns)/log(10),10))); 
cmap = winter(length(idx_plot));
figure;
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren','box','off','tickdir','out');
set(gca,'xlim',[tvec(1) tvec(end)]);
plot(tvec,canc_store(idx_plot,:)');
xlabel('time (s)'); ylabel('MG Vm');

% plot cancelled response, negative image, and sensory response
figure;
plot(tvec,A); hold on; plot(tvec,canc_store(end,:)); plot(tvec,neg_i_store(end,:)); 
xlabel('time (s)'); ylabel('MG Vm');
set(gca,'xlim',[tvec(1) tvec(end)],'box','off','tickdir','out');

% degree of cancellation vs iteration
figure;
plot(mse_ratio);
xlabel('iteration'); ylabel('error');
set(gca,'box','off','tickdir','out');

% plot history of some weights
figure;
plot(whist(1:100:end,:));
set(gca,'xscale','log')

max(max(contribs))

