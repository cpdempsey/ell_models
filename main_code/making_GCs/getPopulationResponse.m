%% gather population summary

function [convspikes_full, convspikes_limit, sptimes_firing, ...
    raster_full, summedSpikes, celltypes_firing, max_freq] = ...
    getPopulationResponse(pth,syn_kernel)

% c_w_syn         = @(x) convolve_with_synaptic_kernel(x', syn_kernel.dt,syn_kernel.tau_on, syn_kernel.tau_off);
% c_w_syn_per         = @(x) convolve_with_synaptic_kernel_pad_sides(x', syn_kernel.dt,syn_kernel.tau_on, syn_kernel.tau_off);
c_w_syn         = @(x) convolve_with_synaptic_kernel_pad_sides(x', syn_kernel.dt,syn_kernel.tau_on, syn_kernel.tau_off);

d=dir([pth '/*GCs*']);

celltypes_firing = [];
summedSpikes = [];
sptimes_firing = [];

waitbar(0);

for ii=1:length(d)
    
    load([pth '/' d(ii).name]);

    dt = GC_models_batch{1}.dt;

    idx_spike = any(meanSpikes');
    current_celltypes = cellnames(idx_spike);
    
    is_periodic = cell2mat(cellfun(@(x) ~isempty(strfind(x,'pause') | strfind(x,'tonic')),...
        cellnames,'uniformoutput',false));
    
    celltypes_firing        =   [celltypes_firing ; current_celltypes];
    
    %     c_w_syn         = @(x) convolve_with_synaptic_kernel_pad_start(x', syn_kernel.dt,...
    %         syn_kernel.tau_on, syn_kernel.tau_off,'padding_constant',);
    
    %     for jj=1:size(meanSpike(idx_spike,:),1)
    %         convspikes(jj,:) = c_w_syn(full(meanSpikes(idx_spike(jj),:)));
%     %     end
%     convspikes = zeros(size(meanSpikes,2),nnz(idx_spike));
%     convspikes(:,is_periodic(idx_spike)) = c_w_syn_per(full(meanSpikes(is_periodic(idx_spike),:)));
%     convspikes(:,~is_periodic(idx_spike)) = c_w_syn(full(meanSpikes(~is_periodic(idx_spike),:)));
%     
    convspikes = c_w_syn(full(meanSpikes(idx_spike,:)));
    
    if ii==1
        convspikes_full = convspikes;
        sptimes_firing = sptimes(idx_spike);
        raster_full = meanSpikes(idx_spike,:);
        summedSpikes = sum(meanSpikes(idx_spike,:),2);
    else
        convspikes_full = [convspikes_full convspikes];
        sptimes_firing = [sptimes_firing ; sptimes(idx_spike)];
        raster_full = [raster_full ; meanSpikes(idx_spike,:)];
        summedSpikes = [ summedSpikes; sum(meanSpikes(idx_spike,:),2)];
    end
    
    waitbar(ii/length(d));
    
end

% frequency limited basis to speed up simulation of learning
[convspikes_limit, max_freq] = freqLimit(convspikes_full,dt/1000);

% save([pth '/population_summary.mat'],'convspikes_full','convspikes_limit','sptimes_firing','raster_full','celltypes_firing','summedSpikes','max_freq','-v7.3');
