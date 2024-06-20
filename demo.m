%% prepare training data and test data
clc
clear
close all
addpath condition_driven_mode_division/
addpath base_model/
addpath utils/
load('data/d11_te.mat')
data = data';
data_train = data(:,[9,10,51]); % samples*variables
data_train_diff = diff(data_train);
data_train(1,:)=[];
data_test = data(:,[9,10,51]);
data_test_diff = diff(data_test);
data_test(1,:)=[];
%% parameter setting 
% Important Note: Please change these parameters according to your needs. 
model_directory_and_name = 'trained_model/divided_models.mat';%the directory to
                       %save the model(including the name of the model)
confidence_level = 0.99; %confidence level of the control limits
alpha = 1.1; % Relaxing factor. The larger the alpha is, the less number of mode will be divided.
lower_limit = 0; % the samples with the value of indicate_varibale lower 
%                 than lower_limit will be discarded, default 0.
indicate_variable = 3;% indicates which column of the data matrix is the
%                       indicate variable, default value is 3, i.e.,...
%                       divide based on the third variable.(指示变量标签)

%% initialize model
% Important Note: You can change base model according to your needs. 
% If so, you should prepare your own class based on the 'Base_model/SFA_class.m' .

% SFA as base model
thresholdv = 1e-3; % 白化阈值 the whitening threshold
monitoring_statistic_num = 4;% the number of monitoring statistic
base_model = SFA_class(thresholdv,confidence_level, monitoring_statistic_num);

% Uncomment the following code if you want to use PCA as base model
% %PCA as base model
% thresholdv = 95; % 白化阈值 the whitening threshold
% monitoring_statistic_num = 2;% the number of monitoring statistic
% base_model = PCA_class(thresholdv,confidence_level, monitoring_statistic_num);
%% modeling
divided_models = condition_driven_mode_division_modeling(data_train, data_train_diff, lower_limit, indicate_variable,...
                alpha, base_model);% the output is a struct named as 'divided_models'
save(model_directory_and_name, 'divided_models'); 
%% monitoring
% Important Note: This is an example showing how to use divided models for monitoring.
% You can use divided models for other purpose. 
% If so, you should replace the 'utils/monitoring.m' function with your own function. 
load(model_directory_and_name, 'divided_models');
base_model.feature_nums = divided_models.feature_nums;
[monitoring_statistics, ctrl_limits, divided_models] = monitoring(data_train, data_train_diff, data_test,data_test_diff, divided_models, base_model);
save(model_directory_and_name, 'divided_models'); 
%% show monitoring results
figure()
plot(data_train(:,indicate_variable)),hold on
for i = 1:size(divided_models.border_segment,1)
    plot([0, size(data_train,1)],[divided_models.border_segment(i),divided_models.border_segment(i)],'--r'), hold on
end
figure()
for i = 1:monitoring_statistic_num
    subplot(monitoring_statistic_num,1,i)
    plot(monitoring_statistics(:,i),'-b'), hold on
    plot(ctrl_limits(:,i),'--k')
end