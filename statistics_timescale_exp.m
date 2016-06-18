function [sta] = statistics_timescale_exp(para_est, aggr_time)
% homogeneous case. 
% get second order statistics [mean;var;P_nodemand;lag1_autocorr] at
% different time scale.
% exponential-lognormal distribution
% input parameter on second unit. 1*5;
if nargin <2
    
    for i = 1:3600
        if mod(3600,i) ==0
            b(i) = 1;
        end
    end
    aggr_time  = find(b == 1)';
end


n = length(aggr_time);
sta = zeros(4,n);
for i = 1:length(aggr_time)
    t1 = aggr_time(i);
    meanexp_m = t1*(para_est(1)*para_est(3)*para_est(5));
    var_exp_m =  (t1./para_est(3)-1+exp(-t1./para_est(3))) .* (2*para_est(5).*(para_est(3)).^3.*(para_est(1).^2+para_est(2).^2));
    Pexp_m = exp(-(para_est(3)+t1)*para_est(5));
    autocorr_m = 0.5*(1-exp(-1/para_est(3)*t1))^2/(1/para_est(3)*t1-1+exp(-1/para_est(3)*t1));
    sta(:,i) = [meanexp_m;var_exp_m;Pexp_m;autocorr_m];
end