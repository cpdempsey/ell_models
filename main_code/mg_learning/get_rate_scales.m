function w = get_rate_scales(convolved_rates)

popmean = mean(abs(convolved_rates),2);
ncells = size(convolved_rates,2);
for i=1:ncells
    [~,ind] = max(abs(convolved_rates(:,i)));
    w(i) = 1/(popmean(ind));
end

w=w/std(w);