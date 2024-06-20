classdef PCA_class
    % principal component analysis

    properties
         pc_num % the number of principal components 
         feature_nums %  equals to the number of principal components
         thresholdv % the whitening threshold
         confidence_level % confidence level of control limit
         monitoring_statistic_num % the number of monitoring statistic
    end

    methods
        function obj = PCA_class(thresholdv,confidence_level, monitoring_statistic_num)
         if nargin > 0
            obj.thresholdv = thresholdv;
            obj.confidence_level = confidence_level;
            obj.monitoring_statistic_num = monitoring_statistic_num;
         end
        end

        function [P, T, S, cov_features] = train(obj, normalized_data, ~)
        % training function for PCA
        % inputs:
        %        normalized_data(samples*variables): data matrix having zero mean
        % outputs:
        %        P:the coefficients from normalized_data to T
        %        T(samples*variables): principal components
        %        S: reconstruct error
        %        cov_features: a cell including covariance of T
            [coeff,score,latent,tsquared,explained] = pca(normalized_data);
            if obj.thresholdv < 1
                pcsnum = sum(latent.^0.5 >= obj.thresholdv);%thresholdv is the threshold value of standard deviation
            else
                pcsnum = sum(cumsum(explained)<obj.thresholdv) + 1;
            end
            P = coeff(:,1:pcsnum);
            T = normalized_data*P;
            S = normalized_data-T*P';
            cov_features = cell(1,2);
            cov_features{1,1} = cov(T);
        end

        function statistics...
            = cal_monitoring_statistics(obj, normalized_data, data_diff, P, cov_features)
            % calculate monitoring statistics
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         data_diff(samples*variables)： the corresponding firt order difference term for
            %               each sample.
            %         P: the coefficients from normalized_data to T
            %         cov_features: covariance of T
            % output: statistics： monitoring statistics.
            T = normalized_data*P;
            S = normalized_data - T*P';
            n = size(T, 1);
            T2 = zeros(n, 1);
            SPE = zeros(n, 1);
            for i = 1:n
                T2(i) = T(i, :)*cov_features{1,1}*T(i, :)';
                SPE(i) = S(i, :)*S(i, :)';
            end
            statistics = [T2, SPE];
        end

        function ctrl_limit = cal_dividing_statistics(obj, normalized_data, P)
            % calculate the control limit of dividing statistic
            % input:
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         P：projection matrix.
            % output:
            %         ctrl_limit： the corresponding control limit of dividing statistic.
            % Important Note: you should choose statistic that monitoring the
            %       relationship between variables as the dividing
            %       statistic.
            T = normalized_data*P(1:obj.pc_num);
            error = normalized_data - T*P(1:obj.pc_num)';
            n = size(error, 1);
            SPE = zeros(n, 1);
            for i = 1:n
                SPE(i) = SPE(i, :)*SPE(i, :)';
            end
            ctrl_limit = ctrl_limit_compute(SPE, obj.confidence_level);
        end

        function [pc_num, sf_num] = find_pcs_num(obj, normalized_data, data_diff, T, ~)
            % find the number of principal components for data slice
            % inputs: 
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         data_diff(samples*variables)： the corresponding firt order difference term for
            %               each sample.
            %         T：latant features;
            % outputs:
            %         pc_num: the number of principal components retained.
            %         sf_num: equals to pc_num.
            pc_num = size(T, 2);
            sf_num = pc_num;
        end
    end
end