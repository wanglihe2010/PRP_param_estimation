function [para_est] = para_estimation_g(demand,T)
dvar = var(demand);
dmean = mean(demand);
dauto = autocorr(demand,2);
dcov = dauto(3) * dvar;
dauto = dauto(2);
dprob = length(find(demand==0))/length(demand);
options = optimset('Algorithm','interior-point','Hessian','bfgs','Display','off');
l = gleastsquare([1,1,1,0.1], dmean, dvar, dauto,dcov, T);
[x,fval,existflag] = fmincon('gleastsquare',[ 0.0950    0.0704   36.0647    0.0062],[],[],[],[],[0 0 0 0],[],[],options,dmean,dvar,dauto,dcov,T);
%[x,fval,existflag] = fmincon('gleastsquare',[ 0.0950    0.0704   36.0647    0.0062],[],[],[],[],[0 0 0 0],[],[],options,dmean,dvar,dauto,dprob,T);

para_est = [x(1:3), x(3), x(4)];
