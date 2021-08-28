%% prepare training data and test data
clc
clear
close all
load('data.mat') % the data set should includes a normal and an abnormal data matrix
data_train = normal;
data_train_diff = diff(data_train);
data_train(1,:)=[];
data_test = abnormal;
%% modeling and monitoring
indicate_variable = 1; %indicates which column of the data matrix the
                          %indicate variable is located.Please change this parameters
                          %according to your needs. 
qieruzhi = 20; %the samples with the value of indicate_varibale lower 
                 %than qieruzhi will be discarded.Please change this parameters
                 %according to your needs. 
confidence_level = 0.99; %confidence level of the control limits
alpha = 1.5; %relaxing factor.Please change this parameters according to your needs. 
model_directory_and_name = 'C:\Users\Desktop\model.mat';%the directory to
                       %save the model(including the name of the model)
                       %Please change this parameters according to your needs. 
main_modeling(data_train,data_train_diff,qieruzhi,indicate_variable,confidence_level,alpha,model_directory_and_name)
[BID, BID_combined, BID_ctrl_limit, BID_ctrl_limit_combined]=main_monitoring(data_test, model_directory_and_name);