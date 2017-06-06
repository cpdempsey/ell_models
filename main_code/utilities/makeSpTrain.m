function sptrain = makeSpTrain(sptimes,time_vec)

dt = time_vec(2) - time_vec(1);

if iscell(sptimes)
    sptrain = cellfun(@(x) histc(x,time_vec)/dt,sptimes,'uniformoutput',false);
    [temp{1:nnz(cellfun(@isempty,sptrain))}]=deal(zeros(1,length(time_vec)));
    sptrain(cellfun(@isempty,sptrain)) = temp;
else
    sptrain = histc(sptimes,time_vec);
    sptrain = sptrain / dt;
    if isempty(sptimes)
        sptrain =  zeros(length(time_vec),1);
    end
end

end

