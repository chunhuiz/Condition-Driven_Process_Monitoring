
function divided_models= condition_driven_mode_division_modeling(data_train, data_train_diff, lower_limit,indicate_variable,...
                alpha, base_model)
% Condition driven mode division modeling function
% inputs:
%       data_train(samples*variables): training data
%       data_train_diff(samples*variables): first order difference of
%       training data
%       lower_limit: the samples with the value of indicate_varibale lower 
%                 than lower_limit will be discarded.
%       indicate_variable: indicates which column of the data matrix the
%                          indicate variable is located.(指示变量标签)
%       alpha: relaxing factor.(宽松系数 )
% outputs: divided_models: a structure including following parameters
%       divided_models.border_segment: the borders of condition segments（每一段的边界）
%       divided_models.border_slice: the borders of condition slices（每一片的边界）
%       divided_models.feature_nums: the number of features（特征个数）
%       divided_models.data_mean_cell: the mean vectors of condition slices（每一片的均值）
%% 判断输入是否合法 assign confidence_level automatically
% if ~isvector(x)
%     error('Input must be a vector')
% end
if ~exist('confidence_level', 'var') || isempty(confidence_level)
    confidence_level = 0.99;
end

%% 数据预处理 data preprocessing
% data_train = fillmissing(data_train,'linear');
% data_train_diff = diff(data_train);
% data_train(1,:)=[];

[num_samples, num_variables] = size(data_train);
% 超参数设定 set the values of hyper parameters
length_min = (num_variables-1)*3*5;% 设置基础数据单元最短长度 set the minimum
                              % number of samples in each condition slices
N = floor(num_samples/(0.01*length_min));
% 去除指示变量值低于切入值的样本  discard the samples with the value of
% indicate_varibale lower than lower_limit
power = data_train(:, indicate_variable);  % 设置划分变量。
id = find(power<lower_limit);
data_train(id,:) = [];
data_train_diff(id,:) = [];
power = data_train(:, indicate_variable);

%% 划分小片 divide original data matrix into different data slices
power_max = max(power);
power_min = min(power);
delta_power = (power_max - power_min)/N;
data_slices = cell(N,1);         % the resulting data slices are stored in a cell
data_diff_slices = cell(N,1);
indicate_variable_slices = cell(N,1);
length = zeros(N,1);
for i = 1:N
    k = find(power>=(power_min+(i-1)*delta_power) & power<(power_min+i*delta_power));
    length(i) = size(k, 1);
    data_slices{i} = data_train(k,[1:indicate_variable-1, indicate_variable+1:end]);
                  %the indicate variable will not be considered when modeling
    data_diff_slices{i} = data_train_diff(k, [1:indicate_variable-1, indicate_variable+1:end]);
    indicate_variable_slices{i} = data_train(k,indicate_variable);
end

%% data slices combinations
% combine adjacent data slices such that the number of samples in each data
% slices is larger than the minimum requirement.
i = 1;
j = 1;
m = 1;
while j < N
    
    if i >= size(data_slices,1)
        break
    end
    
    if size(data_slices{i}, 1)<length_min
        data_slices{i} = [data_slices{i};data_slices{i+1}];
        data_slices(i+1) = [];
        data_diff_slices{i} = [data_diff_slices{i};data_diff_slices{i+1}];
        data_diff_slices(i+1) = [];
        indicate_variable_slices{i} = [indicate_variable_slices{i};indicate_variable_slices{i+1}];
        indicate_variable_slices(i+1) = [];
        j = j + 1;
        m = m + 1;
    else
        i = i + 1;
        j = j + 1;
    end
end
% 如果最后一个基础数据单元矩阵样本数量小于最小样本数，则将这个基础数据矩阵与前一个合并。
if size(data_slices{end},1)<length_min
    data_slices{end-1} = [data_slices{end-1};data_slices{end}];
    data_slices(end) = [];
    data_diff_slices{end-1} = [data_diff_slices{end-1};data_diff_slices{end}];
    data_diff_slices(end) = [];
    indicate_variable_slices{end-1} = [indicate_variable_slices{end-1};indicate_variable_slices{end}];
    indicate_variable_slices(end) = [];
end

%% 基础数据单元分析 analyze each data slice
% 将每个基础单元数据矩阵标准化
N = size(data_slices, 1); %调整后的基础数据单元个数。the number of data slices
                          %after conbination
data_slices_normalized = cell(N,1); %data_diff_slices_normalized = cell(N,1);
data_mean_cell = cell(N,1);  %data_std_cell = cell(N,1);
pc_nums = zeros(N, 1);
sf_nums = zeros(N, 1);
border_slice = zeros(N+1,1);
border_slice(1) = power_min;
for i = 1:N
    border_slice(i+1) = max(indicate_variable_slices{i});
    data_mean_cell{i} = mean(data_slices{i});
    %    data_std_cell{i} = std(data_slices{i});
    data_slices_normalized{i} = (data_slices{i}-repmat(data_mean_cell{i},size(data_slices{i},1),1));%./data_std_cell{i};
    %    data_diff_slices_normalized{i} = data_diff_slices{i}./data_std_cell{i};
    [~, T, S] = base_model.train(data_slices_normalized{i}, data_diff_slices{i});
    [pc_nums(i), sf_nums(i)] = base_model.find_pcs_num(data_slices_normalized{i}, data_diff_slices{i}, T, S);
end

%% 找到出现次数最多的主元个数和慢特征个数
% Find the number of principal components and the number of slow features that 
% appear the most frequently  
table1 = tabulate(pc_nums);
maxcount = max(table1(:,2));
[row1,~]=find(table1(:,2)==maxcount);
pc_num = table1(row1,1);
pc_num = pc_num(1);

table2 = tabulate(sf_nums);
maxcount = max(table2(:,2));
[row2,~]=find(table2(:, 2)==maxcount);
feature_nums = table2(row2,1);
feature_nums = feature_nums(1);

base_model.pc_num = pc_num;
base_model.feature_nums = feature_nums;
%% 计算条件片控制限 estimate control limit of dividing statistic for each condition slices
ctrl_T2s = zeros(N,1);
for i = 1:N
    W = base_model.train(data_slices_normalized{i}, data_diff_slices{i});
    ctrl_T2s(i) = base_model.cal_dividing_statistics(data_slices_normalized{i}, W);
end
%% 子阶段划分 divide condition slices into different condition segments
N  = size(data_slices, 1);
data_segments = data_slices;
data_segments_normalized = data_slices_normalized;
data_diff_segments = data_diff_slices;
indicate_variable_segments = indicate_variable_slices;
label = 1; %标签 segment label
i = 1;
k = 1;

while i <= N
    j = 1;
    while j <= N-i
        % determine whether the control limit of the condition segment model is
        % larger than that of the condition slice model multiplied by alpha.
        L = 3;
        data_unfold = [];
        data_unfold_diff = [];
        for l = i:i+j
            data_unfold = [data_unfold; data_slices_normalized{l}];
            data_unfold_diff = [data_unfold_diff; data_diff_slices{l}];
        end
        W = base_model.train(data_unfold, data_unfold_diff);
        % 比较控制限 compare control limits
        ctrl_T2s_vk = zeros(j+1,1);
        for l = i:i+j
            ctrl_T2s_vk(l-i+1) = base_model.cal_dividing_statistics(data_slices_normalized{l}, W);
        end
        alphas1 = (ctrl_T2s_vk./ctrl_T2s(i:i+j)) >= alpha;
        % if alphas1 >= alpha in 3 consecutive slices, a new segment should
        % be created.
        % 在连续三个小片上大于alpha即为不能合并
        M1 = movsum(alphas1,[L-1 0]);
        if max(M1)>=3
            result = 1;
        else
            result = 0;
        end

        if result == 1       % result==1 meaning that the i+j th slice cannot 
            % be divided into current segment 意味着第i+j个slice不能划入当前子阶段
            label = label + 1; % set new label 设置新的label值，为第i+j个slice打上新label
            k = k + 1;
            break % start a new loop from i+j th slice 从第i+j个slice开始新的循环。
        else      %result==0 meaning that the i+j th slice should be divided
            % into current segment 意味着第i+j个slice是处在当前子阶段.
            % set the same label, merge into current segment 打上相同标签、并入当前子阶段。
%             data_segments_normalized{k+1}(:, num_variables+1) = label;
            data_segments_normalized{k} = [data_segments_normalized{k};data_segments_normalized{k+1}];
            data_segments_normalized(k+1) = [];
%             data_segments{k+1}(:, num_variables+1) = label;
            data_segments{k} = [data_segments{k};data_segments{k+1}];
            data_segments(k+1) = [];
%             data_diff_segments{k+1}(:, num_variables+1) = label;
            data_diff_segments{k} = [data_diff_segments{k};data_diff_segments{k+1}];
            data_diff_segments(k+1) = [];
            indicate_variable_segments{k} = [indicate_variable_segments{k};indicate_variable_segments{k+1}];
            indicate_variable_segments(k+1) = [];
            j = j + 1;
        end
    end
    i = i + j
end

%% 确定划分边界 determine the borders of segments
n = size(indicate_variable_segments, 1);
border_segment = zeros(n+1,1);
border_segment(1) = power_min;
for i = 1:n
    border_segment(i+1) = max(indicate_variable_segments{i});
end

%% save parameters
divided_models = struct;
divided_models.border_segment = border_segment;
divided_models.border_slice = border_slice;
divided_models.indicate_variable = indicate_variable;
divided_models.feature_nums = feature_nums;
divided_models.data_mean_cell = data_mean_cell;
end


