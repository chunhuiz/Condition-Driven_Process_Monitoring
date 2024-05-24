function [BID, BID_combined, BID_ctrl_limit, BID_ctrl_limit_combined]=main_monitoring(data_test, model_directory_and_name)
%% ÔÚÏß¼à²â online monitoring
%   inputs£º
%       data_test(samples*variables)£º test data
%       model_directory_and_name: the directory to save the model(including
%                   the name of tha model) e.g. 'C:\Users\Desktop\model.mat' 
%   outputs:
%       BID(samples-1*4): BID monitoring indices
%       BID_combined(samples-1*1): combined BID monitoring index
%       BID_ctrl_limit(samples-1*4): corresponding control limit for each
%                                    samples and each monitoring indices.
%       BID_ctrl_limit_combined(samples-1*1):corresponding control limit
%                                           for BID_combined.

data_test_diff = diff(data_test);
data_test(1,:) = [];
load(model_directory_and_name, 'model');
[BID, BID_combined, BID_ctrl_limit, BID_ctrl_limit_combined, ~, ~] = monitoring(data_test,data_test_diff,model.border, model.border_slice, model.sf_nums,...
    model.indicate_variable, model.data_mean_cell, model.W_cell, model.GMMmodel_cell,...
    model.BID_ctrl_limit_matrix, model.BID_ctrl_limit_matrix_combined, model.Cont_BID_CL_cell);

end