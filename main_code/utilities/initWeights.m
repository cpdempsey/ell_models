function  wd_init = initWeights(drop_limit,ubs,lbs,lambda,target)

options = optimoptions('fmincon','GradObj','on','display','iter','maxiter',50);

% drop
G = drop_limit;
w0 = 1+randn(size(G,2),1);
c = target*ones(size(G,1),1);
ub = ubs*ones(size(w0));
lb = lbs*ones(size(w0));

GtG = G'*G;
Gtc = G'*c;
ctc = c'*c;
wd_init = fmincon(@(w) L2wAbsTotalVariation(GtG,Gtc,ctc,w,lambda),w0,[],[],[],[],lb,ub,[],options);

end