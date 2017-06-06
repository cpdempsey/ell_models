function [amp_response, amp_noise, dt] = get_amp_response(max_freq)

load 'amp_cell_data.mat';
amp_timebins = amp_timebins/1000;

dt = (amp_timebins(2)-amp_timebins(1));

tran_extend = -.025+dt:dt:.2;

padfront = find(tran_extend>amp_timebins(1),1);
padend = length(tran_extend)-find(tran_extend>amp_timebins(end),1);
tEOD = find(tran_extend>0,1);

amp_response = [amp_response(1)*ones(1,padfront) amp_response amp_response(1)*ones(1,padend)] * dt;
amp_response = amp_response - mean(amp_response(1:tEOD-10));


%convolve with MG cell synaptic kernel!!!
MG_cell.dt = 0.05;
MG_cell.tau_m = 5;
MG_cell.tau_s = 10;

Ps   = convolve_matrix_by_tau(MG_cell.dt,MG_cell.tau_s,amp_response);
amp_response = convolve_matrix_by_tau(MG_cell.dt,MG_cell.tau_m,Ps);
    
amp_response = amp_response/max(abs(amp_response));

amp_response = fft(amp_response);

freq = linspace(0, 1, length(amp_response)+1) / dt;
freq = freq(1:end-1);
freq(freq > 0.5 / dt) = freq(freq > 0.5 / dt) - 1/dt;

amp_response = amp_response(abs(freq) < max_freq);



tau_on = MG_cell.tau_m / 1000;
tau_off = MG_cell.tau_s / 1000;
tran = tran_extend - min(tran_extend);

epsp_shape = 1/(tau_off - tau_on) * (exp(-tran/tau_off) - exp(-tran/tau_on));
epsp_shape = epsp_shape / max(abs(epsp_shape));
ft_epsp = fft(epsp_shape);

freq = linspace(0,1,length(epsp_shape)+1) / dt;
freq = freq(1:end-1);
freq(freq>0.5 / dt) = freq(freq>0.5 / dt) - 1/dt;
ft_epsp = ft_epsp(abs(freq)<max_freq);

amp_noise=zeros(size(ft_epsp));
% freq = 5;
amp_noise(1:end)=1;
% amp_noise(1)=0;
amp_noise = amp_noise.*abs(ft_epsp); %smooth noise by synaptic kernel!
amp_noise = amp_noise * sum(abs(amp_response))/sum(abs(amp_noise));

phase=rand(floor(length(amp_noise)/2),1)*2*pi; %give components random phase
phr=[0; phase; -flipud(phase)]';

amp_noise = amp_noise.*(cos(phr) + 1i*sin(phr));

amp_response = ifft(amp_response);

end
