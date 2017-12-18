function [convolved_rates , synaptic_kernel] = convolve_with_synaptic_kernel_pad_sides(input_activity, dt, tau_on, tau_off,varargin)

% Construct synaptic kernel
tau_slow = max(tau_on, tau_off);
filter_times = (-10*tau_slow):dt:(10*tau_slow);
synaptic_kernel = (filter_times > 0) .* (1 - exp(-filter_times/tau_on)) .* exp(-filter_times/tau_off);
synaptic_kernel = synaptic_kernel / max(synaptic_kernel);

% convolve firing rates with synaptic kernel
convolved_rates = zeros(size(input_activity));
for i=1:size(input_activity,2)
    
    temp=[ones(4000,1)*mean(input_activity(370:460,i)); input_activity(:,i); ones(4000,1)*mean(input_activity(end-100:end,i))];
    temp=conv(temp,synaptic_kernel,'same');
    convolved_rates(:,i)=temp(4001:8500);
    
end

end
