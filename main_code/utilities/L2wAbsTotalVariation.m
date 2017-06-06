% evaluate ||Gw - c||_2^2 + lambda*TV(w)
% where TV(x) is the total variation of x

function [L,grad] = L2wAbsTotalVariation(GtG,Gtc,ctc,w,lambda)

    L = w'*GtG*w - 2*w'*Gtc + ctc + lambda*totalVariation(abs(w));
    
    grad = GtG*w - Gtc + lambda*(length(w)*w - sum(abs(w)).*sign(w));

end

function TV = totalVariation(w)

TV = 0;

for ii=1:length(w)
    for jj=1:length(w)
        
        TV = TV + (w(ii) - w(jj))^2;
        
    end
end



end