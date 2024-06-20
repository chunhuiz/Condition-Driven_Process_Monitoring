function [monitoring_statistics, ctrl_limit, divided_models] = monitoring(data_train, data_train_diff, data_test,data_test_diff,divided_models, base_model)

% inputs: 
%   data_train(training data)
%   data_train_diff(first order difference of training data)
%   data_test(test data)
%   data_test_diff(first order difference of test data)
% outputs：
%   monitoring_statistics（monitoring statistics)
%   ctrl_limit（control limits for monitoring statistics)
%   divided_models (models for each segment)
% *please refer to the comments of modeling if you want to know the meaning
% of each parameters.

%% analyze each segment of training data

num_slice = size(divided_models.border_slice,1)-1;
num_segments = size(divided_models.border_segment,1)-1;
power = data_train(:,divided_models.indicate_variable);
for i = 1:num_slice
    k = find(power>=divided_models.border_slice(i) & power<divided_models.border_slice(i+1));
    data_train(k,[1:divided_models.indicate_variable-1, divided_models.indicate_variable+1:end])...
        = data_train(k,[1:divided_models.indicate_variable-1, divided_models.indicate_variable+1:end]) - divided_models.data_mean_cell{i};
                  %the indicate variable will not be considered when modeling
end

rng(42) %复现结果 for reproducing the results
W_cell = cell(num_segments,1);
cov_features_cell = cell(num_segments,1);
ctrl_limit_matrix = zeros(num_segments,base_model.monitoring_statistic_num);
for i = 1:num_segments
    k = find(power>=divided_models.border_segment(i) & power<divided_models.border_segment(i+1));
    X_train = data_train(k,[1:divided_models.indicate_variable-1, divided_models.indicate_variable+1:end]);
                  %the indicate variable will not be considered when modeling
    X_train_diff = data_train_diff(k, [1:divided_models.indicate_variable-1, divided_models.indicate_variable+1:end]);
    [W, ~, ~, cov_features] = base_model.train(X_train, X_train_diff);
    W_cell{i} = W;
    cov_features_cell{i} = cov_features;
    statistics...
            = base_model.cal_monitoring_statistics(X_train, X_train_diff, W, cov_features);
    ctrl_limits = ctrl_limit_compute(statistics, base_model.confidence_level);
    ctrl_limit_matrix(i, :) = ctrl_limits;
end

divided_models.W_cell = W_cell;
divided_models.cov_features_cell = cov_features_cell;
divided_models.ctrl_limit_matrix = ctrl_limit_matrix;
%       divided_models.W_cell: the coefficient matrices of condition segments（SFA变换矩阵）
%       divided_models.cov_features_cell: the covariance matrix used to
%                                           calculate monitoring statistics
%       divided_models.ctrl_limit_matrix: the control limits for each segment.

%% online monitoring of test data

X_test = data_test;
X_test_diff = data_test_diff;
indicate_variable = divided_models.indicate_variable;
border_segment = divided_models.border_segment;
cov_features_cell = divided_models.cov_features_cell;
[n, ~] = size(X_test);

monitoring_statistics = zeros(n,base_model.monitoring_statistic_num); 
ctrl_limit = zeros(n, base_model.monitoring_statistic_num);

for i=1:n
    if X_test(i, indicate_variable)>= max(border_segment) || X_test(i, indicate_variable)< min(border_segment)
        monitoring_statistics(i,:) = -ones(1,base_model.monitoring_statistic_num);
        ctrl_limit(i,:) = zeros(1,base_model.monitoring_statistic_num);
    else
        
        J = border_segment - X_test(i, indicate_variable);%determine which condition
                        % segments the sample belong to so that corresponding model is called. 确定样本对应哪个功率段以调用模型
        j = sum(J <= 0);
        Y =  divided_models.border_slice- X_test(i, indicate_variable); %%determine which condition
                        % slice the sample belong to so that corresponding mean vector is used to standardize the sample. 确定样本对应哪个功率片以进行标准化
        y = sum(Y <= 0);
        
        X_test_scaled = (X_test(i,[1:indicate_variable-1, indicate_variable+1:end])-divided_models.data_mean_cell{y});
        X_test_scaled_diff = X_test_diff(i, [1:indicate_variable-1, indicate_variable+1:end]);
        cov_features = cov_features_cell{j};
        statistics = base_model.cal_monitoring_statistics(X_test_scaled, X_test_scaled_diff, divided_models.W_cell{j}, cov_features);
        monitoring_statistics(i,:) = statistics;
        ctrl_limit(i, :) = divided_models.ctrl_limit_matrix(j, :);
    end
end

end