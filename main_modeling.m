
function main_modeling(data_train,data_train_diff,qieruzhi,indicate_variable,confidence_level, alpha, model_directory_and_name)
%% 建模 modeling
% inputs:
%       data_train(samples*variables): training data
%       data_train_diff(samples*variables): first order difference of
%                                           training data
%       qieruzhi: the samples with the value of indicate_varibale lower 
%                 than qieruzhi will be discarded.
%       indicate_variable: indicates which column of the data matrix the
%                          indicate variable is located.(指示变量标签)
%       confidence_level: confidence level of control limit.（控制限置信度）
%       alpha: relaxing factor.
%       model_directory_and_name: the directory to save the model(including
%                   the name of the model)e.g. 'C:\Users\Desktop\model.mat'

% the output is a file named as 'model.mat'

[border, border_slice, sf_nums,data_mean_cell, W_cell, GMMmodel_cell,...
             BID_ctrl_limit_matrix, BID_ctrl_limit_matrix_combined, Cont_BID_CL_cell]...
                = modeling(data_train, data_train_diff, qieruzhi, indicate_variable, confidence_level, alpha);
model = struct;
model.border = border;
model.border_slice = border_slice;
model.indicate_variable = indicate_variable;
model.sf_nums = sf_nums;
model.data_mean_cell = data_mean_cell;
model.W_cell = W_cell;
model.GMMmodel_cell = GMMmodel_cell;
model.BID_ctrl_limit_matrix = BID_ctrl_limit_matrix;
model.BID_ctrl_limit_matrix_combined = BID_ctrl_limit_matrix_combined;
model.Cont_BID_CL_cell = Cont_BID_CL_cell;
save(model_directory_and_name, 'model');

end