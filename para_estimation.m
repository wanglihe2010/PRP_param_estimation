function para_est = para_estimation(demand1,demand2,t1,t2)
%optchoice : 1 - get alpha and beta just by equation
%            2 - least square 
optchoice =1;
%--------statistic property -------
var1 = var(demand1)/t1/t1;
var2 = var(demand2)/t2/t2;
mean1 = mean(demand1)/t1;
mean2 = mean(demand2)/t2;
p1 = length(find(demand1==0))/length(demand1);
p2 = length(find(demand2==0))/length(demand2);
if p1 == p2
    a=0;
end

% ---------tau, alpha and lambda estimation ------------
tau_est = (t2*log(p1) - t1 * log(p2))/(log(p2)-log(p1));
lambda_est = log(p1/p2)/(t2-t1);
alpha_est = mean1 / tau_est/lambda_est;
% ----------theta estimation ----------------

theta_est = (t1*var1 - t2 * var2)/(var2-var1);
theta_est = 2*tau_est;

% ---- buchburger's variance -----
beta_est = sqrt(var1*(theta_est + t1)/lambda_est/tau_est/theta_est - alpha_est^2);
% ----- more accurate variance -----
var1 = var(demand1);
var2 = var(demand2);
beta_est = sqrt(var1/(2*lambda_est*(tau_est)^3*(t1/tau_est-1+exp(-t1/tau_est)))- alpha_est^2);


if var1/(2*lambda_est*(tau_est)^3*(t1/tau_est-1+exp(-t1/tau_est)))- alpha_est^2 < 0
    optchoice = 2;
end
% [extension]------------------ least square method -----------------------
 
if optchoice == 2;
    l = leastsquare([0.2, 0.1], var1,var2,mean1,mean2,tau_est,lambda_est,t1,t2);
    if ~isfinite(l)
        c = 0;
    end
    options = optimset('Algorithm','interior-point','Hessian','bfgs','Display','off');
    [x,fval,existflag] = fmincon('leastsquare',[0.2,0.1],[],[],[],[],[0 0 ],[],[],options,var1,var2,mean1,mean2,tau_est,lambda_est,t1,t2);

    alpha_est = x(1);
    beta_est = x(2);
    l = leastsquare(x, var1,var2,mean1,mean2,tau_est,lambda_est,t1,t2)
    
end

para_est = [alpha_est beta_est tau_est theta_est lambda_est];
