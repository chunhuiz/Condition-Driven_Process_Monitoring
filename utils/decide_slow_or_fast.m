%% determine the number of slow features
% inputs: 
%         normalized_data(samples*variables)£º data matrix after normalization.
%         data_diff(N*J)£º the corresponding firt order difference term for
%               each sample.
%         T£ºthe static feature obtained by sfa.
%         S£ºthe corresponding first order difference of T.
% outputs:
%         sf_num£ºthe number of slow features.
function sf_num = decide_slow_or_fast(normalized_data, data_diff, T, S)

    delta_data_criterion = diag(data_diff'*data_diff)./diag(normalized_data'*normalized_data);
    index = floor(size(normalized_data,2)*0.9);
    slow_fast_criterion = sort(delta_data_criterion);
    slow_fast_criterion = slow_fast_criterion(index);
    feature_value = diag(S'*S)./diag(T'*T);
    sf_num = sum(feature_value-slow_fast_criterion<=0);

end
