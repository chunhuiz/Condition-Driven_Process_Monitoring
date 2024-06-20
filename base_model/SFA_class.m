classdef SFA_class
    % slow feature analysis

    properties
         pc_num % the number of principal components 
         feature_nums %  the number of slow features
         thresholdv % the whitening threshold
         confidence_level % confidence level of control limit
         monitoring_statistic_num % the number of monitoring statistic
    end

    methods
        function obj = SFA_class(thresholdv,confidence_level, monitoring_statistic_num)
         if nargin > 0
            obj.thresholdv = thresholdv;
            obj.confidence_level = confidence_level;
            obj.monitoring_statistic_num = monitoring_statistic_num;
         end
        end

        function [W, T, S, cov_features] = train(obj, normalized_data, data_diff)
        % training function for slow feature analysis
        % inputs:
        %        normalized_data(samples*variables): data matrix having zero mean
        %        data_diff(samples*variables): first order difference of normalized_data 
        % outputs:
        %        W: the coefficients from normalized_data to T
        %        T: slow features
        %        S: first order difference of T
        %        cov_features: a 1*2 cell including covariance of Ss and Sf
            [p,l,~]=pcacov(cov(normalized_data));
            if obj.thresholdv < 1
                pcsnum = sum(l.^0.5 >= obj.thresholdv);%thresholdv is the threshold value of standard deviation
            else
                pcsnum = obj.thresholdv;
            end
            
            z=normalized_data*p(:,1:pcsnum)/(diag(l(1:pcsnum).^0.5));   
            Wht=p(:,1:pcsnum)/(diag(l(1:pcsnum).^0.5));
            
            zdiff = data_diff*Wht;
             
            P = pcacov(cov(zdiff));
            P = flip(P, 2);
            T=z*P;
            W=Wht*P;
            S=zdiff*P;

            [~, sf_num] = obj.find_pcs_num(normalized_data, data_diff, T, S);
            cov_features = cell(1,2);
            cov_features{1,1} = cov(S(:,1:sf_num));
            cov_features{1,2} = cov(S(:,sf_num+1:end));
        end

        function statistics...
            = cal_monitoring_statistics(obj, normalized_data, data_diff, W, cov_features)
            % calculate monitoring statistics
            % input:
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         data_diff(samples*variables)： the corresponding firt order difference term for
            %               each sample.
            %         W: the coefficients from normalized_data to T
            %         cov_features: a cell including covariance of Ss and Sf
            % output： 
            %          statistics：monitoring statistics
   
            T = normalized_data*W;
            S = data_diff*W;
            Ts = T(:, 1:obj.feature_nums);
            Tf = T(:, obj.feature_nums+1:end);
            n = size(Ts, 1);
            T2s = zeros(n, 1);
            T2f = zeros(n, 1);
            for i = 1:n
                T2s(i) = Ts(i, :)*Ts(i, :)';
                T2f(i) = Tf(i, :)*Tf(i, :)';
            end
            Ss = S(:, 1:obj.feature_nums);
            Sf = S(:, obj.feature_nums+1:end);
            S2s = diag(Ss * pinv(cov_features{1,1}) * Ss');
            S2f = diag(Sf * pinv(cov_features{1,2}) * Sf');
            statistics = [T2s, T2f, S2s, S2f];
        end

        function ctrl_limit = cal_dividing_statistics(obj, normalized_data, W)
            % calculate dividing statistic% calculate the control limit of dividing statistic
            % input:
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         W：projection matrix.
            % output:
            %         ctrl_limit： the corresponding control limit of dividing statistic.
            % Important Note: you should choose statistic that monitoring the
            %       relationship between variables as the dividing
            %       statistic.
            T = normalized_data*W;
            Ts = T(:, 1:obj.feature_nums);
            n = size(Ts, 1);
            T2s = zeros(n, 1);
            for i = 1:n
                T2s(i) = Ts(i, :)*Ts(i, :)';
            end
            ctrl_limit = ctrl_limit_compute(T2s, obj.confidence_level);
        end

        function [pc_num, sf_num] = find_pcs_num(obj, normalized_data, data_diff, T, S)
            % find the number of principal components and slow features
            % inputs: 
            %         normalized_data(samples*variables)： data matrix after normalization.
            %         data_diff(samples*variables)： the corresponding firt order difference term for
            %               each sample.
            %         T： slow features.
            %         S： the firt order difference of slow features.
            % outputs:
            %         pc_num: the number of principal components retained.
            %         sf_num: the number of slow features.
            pc_num = size(T, 2);
            delta_data_criterion = diag(data_diff'*data_diff)./diag(normalized_data'*normalized_data);
            index = floor(size(normalized_data,2)*0.9);
            slow_fast_criterion = sort(delta_data_criterion);
            slow_fast_criterion = slow_fast_criterion(index);
            feature_value = diag(S'*S)./diag(T'*T);
            sf_num = sum(feature_value-slow_fast_criterion<=0);
        end
    end
end