%% find the number of principal components for data slice
% inputs: 
%         normalized_data(samples*variables)： data matrix after normalization.
%         data_diff(samples*variables)： the corresponding firt order difference term for
%               each sample.
%         thresholdv：a threshold value to determine how many PCs are 
%                    retained for whitening
% outputs:
%         pc_num: the number of principal components retained.
%         sf_num: the number of slow features.
function [pc_num, sf_num] = find_pcs_num(normalized_data, data_diff,thresholdv)

    [T, ~, S, ~] = sfa(normalized_data, data_diff,thresholdv);
    pc_num = size(T, 2);
    %判断快慢特征数量
    sf_num = decide_slow_or_fast(normalized_data, data_diff, T, S);
end
