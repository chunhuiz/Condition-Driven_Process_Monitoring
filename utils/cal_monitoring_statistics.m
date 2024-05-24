%% calculate monitoring indices and conresponding control limits
%         normalized_data(samples*variables)£º data matrix after normalization.
%         data_diff(samples*variables)£º the corresponding firt order difference term for
%               each sample.
%         confidence_level£ºthe confidence level of the control limit eg. 0.99
%         thresholdv£ºa threshold value to determine how many PCs are 
%                    retained for whitening
% Êä³ö²ÎÊý Ts, Tf, Ss, Sf£º static slow feature, static fast feature,
%                       dynamic slow feature, dynamic fast feature.
%          W£ºloading matrix obtained by SFA
%          ctrl_Ts, ctrl_Tf, ctrl_Ss, ctrl_Sf£º the corresponding control
%               limits for four features.
function [Ts, Tf, Ss, Sf, W, ctrl_T2s, ctrl_T2f, ctrl_S2s, ctrl_S2f] ...
    = cal_monitoring_statistics(normalized_data, data_diff, confidence_level, thresholdv, sf_num)
    
    [T, W, S, ~] = sfa(normalized_data, data_diff,thresholdv);
    
    Ts = T(:, 1:sf_num);
    Tf = T(:, sf_num+1:end);
    n = size(Ts, 1);
    T2s = zeros(n, 1);
    T2f = zeros(n, 1);
    for i = 1:n
        T2s(i) = Ts(i, :)*Ts(i, :)';
        T2f(i) = Tf(i, :)*Tf(i, :)';
    end
    ctrl_T2s = ctrl_limit_compute(T2s, confidence_level);
    ctrl_T2f = ctrl_limit_compute(T2f, confidence_level);
    
    Ss = S(:, 1:sf_num);
    Sf = S(:, sf_num+1:end);
    S2s = diag(Ss * pinv(cov(Ss)) * Ss');
    S2f = diag(Sf * pinv(cov(Sf)) * Sf');
    ctrl_S2s = ctrl_limit_compute(S2s, confidence_level);
    ctrl_S2f = ctrl_limit_compute(S2f, confidence_level);
end
