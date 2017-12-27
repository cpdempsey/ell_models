function input = draw_MF_input(mossyraster)

trials_to_choose_from = setdiff(1:size(mossyraster,1),mf_trials);

if isempty(trials_to_choose_from)
    mf_trials_sub = mf_trials;
end

while isempty(trials_to_choose_from)
    [~,mf_idx] = unique(mf_trials_sub);
    mf_trials_sub(mf_idx) = 0;
    trials_to_choose_from = setdiff(1:size(mossyraster,1),mf_trials_sub');
end

trnum = randi(length(trials_to_choose_from));
trnum = trials_to_choose_from(trnum);
input = mossyraster(trnum,:);

% trnum = randi(size(mossyraster,1));
% input = mossyraster(trnum,:);

end

