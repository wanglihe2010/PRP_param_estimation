function para_est = para_estimation(demand1,demand2,t1,t2)
%optchoice : 1 - get alpha and beta just by equation
%            2 - least square 
optchoice =1;
%--------statistic property -------
var1 = var(demand1)/t1/t1;
var2 = var(demand2)/t2/t2;
mean1 = mean(demand1);
mean2 = mean(demand2);
p1 = length(find(demand1==0))/length(demand1);
p2 = length(find(demand2==0))/length(demand2);
if p1 == p2
    a=0;
end

% ---------tau, alpha and lambda estimation ------------
tau_est = (t2*log(p1) - t1 * log(p2))/(log(p2)-log(p1));
lambda_est = log(p1/p2)/(t2-t1);
alpha_est = mean1 / tau_est/lambda_est/t1;
% ----------theta estimation ----------------

theta_est = (t1*var1 - t2 * var2)/(var2-var1);
theta_est = 2*tau_est;

% ---- buchburger's variance -----
beta_est = sqrt(var1*(theta_est + t1)/lambda_est/tau_est/theta_est - alpha_est^2);
% ----- more accurate variance -----
var1 = var(demand1);
var2 = var(demand2);

[alpha_beta] = getalphabeta(mean1, var1, tau_est, lambda_est, t1);


para_est = [alpha_beta tau_est theta_est lambda_est];
