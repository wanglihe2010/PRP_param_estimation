function [para_est] = para_estimation_new(demand,T)
demand_reshape1 = reshape(demand, 24*3600/T, length(demand)/24/3600*T);
demand_reshape= zeros(length(demand)/24,24);
for ii = 1:24
    demand_reshape(:,ii) = reshape(demand_reshape1(3600/T*(ii-1)+1:3600/T*ii,:),length(demand)/24,1);
end
[para_est] = zeros(4,24);
for i = 1:24
    dvar = var(demand_reshape(:,i));
    dmean = mean(demand_reshape(:,i));
    %dauto = autocorr(demand_reshape(:,i),2);
    dauto = cov(demand_reshape(1:end-1,i),demand_reshape(2:end,i))/dvar;
    dauto = dauto(2);
    dprob = length(find(demand_reshape(:,i)==0))/length(demand_reshape(:,i));
    syms x;
    para_est(3,i) = double(solve(dauto - 0.5*(1-exp(-1/x*T))^2/(1/x*T-1+exp(-1/x*T))==0,x));
    para_est(4,i) = -log(dprob)/(para_est(3,i)+T);
    para_est(1,i) = dmean / para_est(3,i)/para_est(4,i)/T;
    %para_est(2,i) = sqrt(dvar/(2*para_est(4,i)*(para_est(3,i))^3*(T/para_est(3,i)-1+exp(-T/para_est(3,i))))- para_est(1,i)^2);
    para_est(1:2,i) = getalphabeta(dmean, dvar, para_est(3,i), para_est(4,i), T)';
    
end



