function result = more_than_alpha(data_slices_normalized, data_diff_slices, i, j, sf_num, ctrl_T2s, ctrl_T2f, alpha, confidence_level, thresholdv)
% determine whether the control limit of the condition segment model is
% larger than that of the condition slice model multiplied by alpha.
N = 3;
data_unfold = [];
data_unfold_diff = [];
for k = i:i+j
    data_unfold = [data_unfold; data_slices_normalized{k}];
    data_unfold_diff = [data_unfold_diff; data_diff_slices{k}];
end
[~, W, ~, ~] = sfa(data_unfold, data_unfold_diff, thresholdv);
%比较控制限
ctrl_T2s_vk = zeros(j+1,1);
ctrl_T2f_vk = zeros(j+1,1);
for k = i:i+j
    T = data_slices_normalized{k}*W;
    Ts = T(:, 1:sf_num);
    Tf = T(:, sf_num+1:end);
    T2s = diag(Ts*Ts');
    T2f = diag(Tf*Tf');
    ctrl_T2s_vk(k-i+1) = ctrl_limit_compute(T2s ,confidence_level);
    ctrl_T2f_vk(k-i+1) = ctrl_limit_compute(T2f ,confidence_level);
end
alphas1 = (ctrl_T2s_vk./ctrl_T2s(i:i+j)) >= alpha;
alphas2 = (ctrl_T2f_vk./ctrl_T2f(i:i+j)) >= alpha;
%连续三个大于即为不能合并
M1 = movsum(alphas1,[N-1 0]);
M2 = movsum(alphas2,[N-1 0]);
if max(M1)>=3 || max(M2)>=3
    result = 1;
else
    result = 0;
end

end
