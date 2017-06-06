function [signals, max_freq,dtn] = freqLimit(signals,dtf)

max_freq = 300;

% frequency limit things for speed:
% dtf = 5e-05;

T = length(signals)*dtf;

freq = linspace(0, 1, size(signals,1)+1) / dtf;
freq = freq(1:end-1);
freq(freq > 0.5 / dtf) = freq(freq > 0.5 / dtf) - 1/dtf;

signals = fft(signals);
signals = signals(abs(freq)<max_freq,:);
signals = ifft(signals);

dtn = T / length(signals);

end