% compute contribute plot for BID using GDC
function Cont_S2 = contribute_plot_GMM_GDC(x_test, x_test_diff, feature,  W_original, sf_nums, k, GMM_model, Pix) 
[x_raw,x_col] = size(x_test);
cluster_num = size(Pix,2);
Cont_S2 = zeros(x_raw,x_col);
inv_W_original = pinv(W_original);
for j = 1:cluster_num
    if k<=2
        Sigma_inv = inv(diag(GMM_model.Sigma(:,:,j)));
        x_orignial = x_test;
    else
        Sigma_inv = inv(GMM_model.Sigma(:,:,j));
        x_orignial = x_test_diff;
    end
    mu = GMM_model.mu(j,:);
    if k == 1 || k ==3
        W = W_original(:,1:sf_nums);
        x_hat = mu*inv_W_original(1:sf_nums,:);
        M = real((W*Sigma_inv*W')^0.5);
    else
        W = W_original(:,sf_nums+1:end);
        x_hat = mu*inv_W_original(sf_nums+1:end,:);
        M = real((W*Sigma_inv*W')^0.5);
    end
    
    Cont_S2 = Pix(:,j).*(((x_orignial-x_hat)*M).^2 )+ Cont_S2;
    
end
% Cont_S2 = Cont_S2*Pix';
% Cont_S2 = Cont_S2.*100./sum(Cont_S2);
end