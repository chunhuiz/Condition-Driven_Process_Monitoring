function [T, W, S, L] = sfa(X, Xdiff,thresholdv)
% slow feature analysis
% inputs:
%        X(samples*variables): data matrix having zero mean
%        Xdiff(samples*variables): first order difference of X 
%        thresholdv: a threshold value to determine how many PCs are 
%                    retained for whitening
%outputs:
%        T(samples*variables): slow features
%        W:the coefficients from X to T
%        S:first order difference of T

[p,l,~]=pcacov(cov(X));
if thresholdv < 1
    pcsnum = sum(l.^0.5 >= thresholdv);%thresholdv is the threshold value of standard deviation
else
    pcsnum = thresholdv;
end

z=X*p(:,1:pcsnum)/(diag(l(1:pcsnum).^0.5));   
Wht=p(:,1:pcsnum)/(diag(l(1:pcsnum).^0.5));

zdiff = Xdiff*Wht;
 
[P,L,~]=pcacov(cov(zdiff));
P = flip(P, 2);
L = flip(L);
T=z*P;
W=Wht*P;
S=zdiff*P;

end
