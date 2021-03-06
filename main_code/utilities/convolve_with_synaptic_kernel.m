function convolved_rates = convolve_with_synaptic_kernel(input_activity, dt, tau_on, tau_off)

% Construct synaptic kernel
tau_slow = max(tau_on, tau_off);
filter_times = (-10*tau_slow):dt:(10*tau_slow);
synaptic_kernel = (filter_times > 0) .* (1 - exp(-filter_times/tau_on)) .* exp(-filter_times/tau_off);
synaptic_kernel = synaptic_kernel / sum(synaptic_kernel);

% convolve firing rates with synaptic kernel
convolved_rates = zeros(size(input_activity));
for i=1:size(input_activity,2)
    
    % suggestion: pad with edge values to keep from adding artifacts to the
    % pause neurons?
    %     temp=[ones(4000,1)*input_activity(1,i); input_activity(:,i)];
    %     temp=conv(temp,synaptic_kernel,'same');
    %     convolved_rates(:,i)=temp(4001:end);
    
    tmpConv = conv([input_activity(:,i); input_activity(:,i)], synaptic_kernel, 'same');
    convolved_rates(:,i) = tmpConv((size(input_activity, 1)+1):end);
end

end
