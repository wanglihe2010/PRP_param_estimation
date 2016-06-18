function l = leastsquare(x, var1_obs,var2_obs,mean1_obs,mean2_obs,tau_est,lambda_est,t1,t2);
alpha = x(1);
beta = x(2);

%var_exp = lambda_est*tau_est*(alpha^2+beta^2)/3600*2*tau_est*60/(2*tau_est*60+t1)*t1*t1;
var_exp1 = 2*lambda_est*(tau_est).^3*(alpha.^2+beta.^2)*(t1/tau_est-1+exp(-t1/tau_est));
var_exp2 = 2*lambda_est*(tau_est).^3*(alpha.^2+beta.^2)*(t2/tau_est-1+exp(-t2/tau_est));

mean_exp1 =alpha * tau_est * lambda_est;
mean_exp2 = alpha * tau_est * lambda_est;
l1 = (mean_exp1 - mean1_obs)^2/mean_exp1^2 +  (var_exp1 - var1_obs)^2/var_exp1^2;
l2 =   (mean_exp2 - mean2_obs)^2/mean_exp2^2  + (var_exp2 - var2_obs)^2/var_exp2^2;
% we can choose consider only one or both T statitics
l = l1+l2;
