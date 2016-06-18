function [time, int, dur, sim_s] = pulse_generation(para,lambda, int_dist, dur_dist)
%%% alpha beta, L/s; tau omega, s; lambda, /s, 1*24; T, s
%%% output: time(s),int(L/s),dur(s), sim_s(s)
clear t demand
%dur_dist = 0; % 1 - lognormal, 0 - exponential, 2 - gamma
%int_dist = 2; % 1 - lognormal, 0 - exponential, 2 - weibull, 3- gamma

% -----------------------set variables and convert units------------------
alpha = para(1,:);
beta = para(2,:);   %convert to L/s
tau = para(3,:);   %convert to second
%omega = para(4,:);

sim_s = 24*2140*3600;
%sim_s = 24*2560*7*10*60;

%-------------- create event happen time line-----------------------------
n = floor(sim_s * mean(lambda));% designed total arrival times
n_new = poissrnd(n);  %acctual total arrival times
lambda_cdf = cumsum(lambda) / sum(lambda); 

choose_day = randi([0 sim_s/3600/24-1],n_new,1);%unit in day
unif_min = rand(1,n_new); %uniform variable for minute selection
%choose_min = getmin(unif_min', lambda_cdf);
lambda_cdf = [0, lambda_cdf]; 
a = bsxfun(@minus,repmat(unif_min',1,25),lambda_cdf);
a(a<0) = 1;
[b, index] = max(-a,[],2);
%index = max(find((p - lambdacdf)>0));
choose_sec = (unif_min' - lambda_cdf(index)')./(lambda_cdf(index+1) - lambda_cdf(index))'*3600 + (index-1)*3600;

time = choose_day*24*3600 + choose_sec;
time = sort(time);
time = round(time);
%---------------------create insensity and duration-----------------------

dur = zeros(n_new,1);
int = zeros(n_new,1);

time_h = floor(time/3600);
for i = 1:24
    index = find(mod(time_h,24) ==i-1);
    % ... for intensity (L/min)
    if int_dist == 0
        int(index) = exprnd(alpha(i),length(index),1);
    elseif int_dist == 1
        mu_a = log(alpha(i)^2/sqrt(beta(i)^2+alpha(i)^2));
        sigma_b = sqrt(log(beta(i)^2/alpha(i)^2+1));
        int(index) = lognrnd(mu_a,sigma_b,length(index),1) ;
    elseif int_dist == 2
        b = (beta(i)/alpha(i))^-1.086;
        a = alpha(i)/gamma(1+1/b);
        int(index) = wblrnd(a,b,length(index),1);
    else
        a = alpha(i)^2/beta(i)^2;
        b = beta(i)^2/alpha(i);
        int(index) = gamrnd(a,b,length(index),1);
    end
    
    if dur_dist == 0
        dur(index) = exprnd(tau(i),length(index),1) ;
    elseif dur_dist == 1
        % ... for duration (min)
        mu_t = log(tau(i)^2/sqrt(omega(i)^2+tau(i)^2));
        sigma_o = sqrt(log(omega(i)^2/tau(i)^2+1));
        dur(index) = lognrnd(mu_t,sigma_o,length(index),1) ;
    else
        a = tau(i)^2/omega(i)^2;
        b = omega(i)^2/tau(i);
        dur(index) = gamrnd(a,b,length(index),1);
    end
end

dur = round(dur);