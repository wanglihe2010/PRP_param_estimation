% non-busy time probability accurancy comparetion on different aggregation
% time scale. 
load('demand_e.mat')

homeint = demand_e{1,2}(:,3);  % unit : L/s
homedur = demand_e{1,2}(:,2);  % unit : s
hometime = demand_e{1,2}(:,1); % unit : s
sim_s = 214*86400;

% --------------Buchberger's parameter ----------

alpha =zeros(1,24);
beta = zeros(1,24);
tau = zeros(1,24);
omega = zeros(1,24);
timehour = mod(hometime,86400);
timehour = floor(timehour/3600);
lambda_24 = hist(timehour,[0:1:23])/214/3600;
for i = 1:24
    index = find(timehour == i-1);
    tau(i) = mean(homedur(index));
    omega(i) = std(homedur(index));
    alpha(i) = sum(homedur(index).*homeint(index))/sum(homedur(index));
    beta(i) = sqrt(sum(homedur(index).*homeint(index).^2)/sum(homedur(index)) - alpha(i) ^2);

    alpha_o(i) = mean(homeint(index));
    beta_o(i) = std(homeint(index));
    
end
para_est_b = [alpha;beta;tau;omega;lambda_24];
homedemand = pulse_aggreation(hometime, homeint, homedur, sim_s, 1);
% convert to 24 hours matrix n*24;
demand_reshape1 = reshape(homedemand, 24*3600/1, length(homedemand)/24/3600*1);
demand_reshape= zeros(214*3600/1,24);
for ii = 1:24
    demand_reshape(:,ii) = reshape(demand_reshape1(3600/1*(ii-1)+1:3600/1*ii,:),214*3600/1,1);
end
for ij = 1:24
    i = 1000;
    a = [];
    while isempty(a)
        autocorrdemand = autocorr(demand_reshape(:,ij),i);
        a = find(autocorrdemand <= 0);
        i = i*10;
    end
    theta(ij) = sum(autocorrdemand(1:a(1)-1));
end




% -------------aggregation time vector ---------------
for i = 1:3600
    if mod(3600,i) ==0
        b(i) = 1;
    end
end
aggr_time  = find(b == 1)';
% ----------------observed second order statitics on different
% aggregation time -------------------
[meanobs, varobs, Pobs, autoobs] = statistics_timescale_obs(hometime, homeint, homedur, sim_s, aggr_time);

%--------------g's parameter estimation ----------------
j=23;  %determine observed data aggregation time
t1 = aggr_time(j);
homedemand1 = pulse_aggreation(hometime, homeint, homedur, sim_s, t1);

demand_reshape1 = reshape(homedemand1, 24*3600/t1, length(homedemand1)/24/3600*t1);
demand_reshape= zeros(214*3600/t1,24);
for ii = 1:24
    demand_reshape(:,ii) = reshape(demand_reshape1(3600/t1*(ii-1)+1:3600/t1*ii,:),214*3600/t1,1);
end
for i = 1:24
[para_est_g(:,i)] = para_estimation_g(demand_reshape(:,i),t1)';

end
% --------------my parameter estimation -------------------
% using t2 with different ratio
n = length(aggr_time);
alpha_24_est = zeros(n,24);
beta_24_est = zeros(n,24);
tau_24_est = zeros(n,24);
lambda_24_est = zeros(n,24);

for i =j+1:n
    t2 = aggr_time(i);
    if mod(t2,t1)==0
        homedemand2 = pulse_aggreation(hometime, homeint, homedur, sim_s, t2);
        [para_est] = para_estimation_24(homedemand1,homedemand2,t1,t2);
        alpha_24_est(i,:) = para_est(:,1)';
        beta_24_est(i,:) = para_est(:,2)';
        tau_24_est(i,:) = para_est(:,3)';
        lambda_24_est(i,:) = para_est(:,5)';
        
    end
end
 
alpha_24_est(1,:) = alpha;
beta_24_est(1,:) = beta;
tau_24_est(1,:) = tau;
lambda_24_est(1,:) = lambda_24;








%-------------expect statistic of Buchberger's parameter -----------
for i = 1:24
    [sta_b] = statistics_timescale_exp(para_est_b(:,i)', aggr_time);
    meanexp_b(:,i) = sta_b(1,:)';
    varexp_b(:,i) = sta_b(2,:)';
    Pexp_b(:,i) = sta_b(3,:)';
    autoexp_b(:,i) = sta_b(4,:)';
    %Pexp(i,:) = exp(-(tau+aggr_time(i)).*lambda_24);
    %meanexp(i,:) = aggr_time(i)*(alpha.*tau.*lambda_24);
    varexp_b_b(:,i) = (lambda_24(i)*(tau(i)*(alpha(i)^2+beta(i)^2)*theta(i)./(theta(i)+aggr_time).*aggr_time.^2));
end

Perror_b = abs(-Pexp_b + Pobs)./Pexp_b*100;
meanerror_b = abs(meanobs - meanexp_b)./meanexp_b *100;
varerror_b = abs(varobs - varexp_b)./varexp_b *100;
varerror_b_b = abs(varobs - varexp_b_b)./varexp_b_b *100;
autoerror_b = abs(autoobs - autoexp_b)./autoexp_b *100;
%-------------expect statistic of g's parameter -----------
for i = 1:24
    [sta_g] = statistics_timescale_exp(para_est_g(:,i)', aggr_time);
    meanexp_g(:,i) = sta_g(1,:)';
    varexp_g(:,i) = sta_g(2,:)';
    Pexp_g(:,i) = sta_g(3,:)';
    autoexp_g(:,i) = sta_g(4,:)';
    %Pexp(i,:) = exp(-(tau+aggr_time(i)).*lambda_24);
    %meanexp(i,:) = aggr_time(i)*(alpha.*tau.*lambda_24);
    %var_exp_b(i,:) = (lambda_24.*(tau.*(alpha.^2+beta.^2)*theta/(theta+aggr_time(i))*aggr_time(i)^2));
end

Perror_g = abs(-Pexp_g + Pobs)./Pexp_g*100;
meanerror_g = abs(meanobs - meanexp_g)./meanexp_g *100;
varerror_g = abs(varobs - varexp_g)./varexp_g *100;
autoerror_g = abs(autoobs - autoexp_g)./autoexp_g *100;


% -------------- compare mean relative abs error of statistics of 24 hours -----------------
p_compare = [[mean(abs(Perror_b),2);1],[mean(abs(Perror_g),2);t1]] ;
m_compare =  [[mean(abs(meanerror_b),2);1],[mean(abs(meanerror_g),2);t1]] ;
v_compare =  [[mean(abs(varerror_b),2);1],[mean(abs(varerror_g),2);t1]] ;
a_compare =  [[mean(abs(autoerror_b),2);1],[mean(abs(autoerror_g),2);t1]] ;
Pexp_temp = zeros(45,24);
meanexp_temp= zeros(45,24);
var_exp_g_temp= zeros(45,24);
autoexp_temp = zeros(45,24);
for i = j+1:n
    if lambda_24_est(i,1) > 0
        
        for ii = 1:n
            Pexp_temp(ii,:) = exp(-(tau_24_est(i,:)+aggr_time(ii)).*lambda_24_est(i,:));
            meanexp_temp(ii,:) = aggr_time(ii)*(alpha_24_est(i,:).*tau_24_est(i,:).*lambda_24_est(i,:));
            var_exp_g_temp(ii,:) =  (aggr_time(ii)./tau_24_est(i,:)-1+exp(-aggr_time(ii)./tau_24_est(i,:))) .* (2*lambda_24_est(i,:).*(tau_24_est(i,:)).^3.*(alpha_24_est(i,:).^2+beta_24_est(i,:).^2));
            autoexp_temp(ii,:) = 0.5*(1-exp(-1./tau_24_est(i,:)*aggr_time(ii))).^2./(1./tau_24_est(i,:)*aggr_time(ii)-1+exp(-1./tau_24_est(i,:)*aggr_time(ii)));
        end
        % expect P 
        Perror_temp = abs(-Pexp_temp + Pobs)./Pexp_temp*100;
        p_compare = [p_compare,[mean(abs(Perror_temp),2);aggr_time(i)]];  
        
        % each coloum is the error for different aggregation time 
        % each row is the error for different T2 in parameter estimation 
        
        % expect mean
        meanerror_temp = abs(meanobs - meanexp_temp)./meanexp_temp *100;
        m_compare = [m_compare, [mean(abs(meanerror_temp),2);aggr_time(i)]];
        % expect var
        varerror_temp = abs(varobs - var_exp_g_temp)./var_exp_g_temp *100;
        v_compare = [v_compare, [mean(abs(varerror_temp),2);aggr_time(i)]];
        % expect auto
        autoerror_temp = abs(autoobs - autoexp_temp)./autoexp_temp *100;
        a_compare = [a_compare, [mean(abs(autoerror_temp),2);aggr_time(i)]];
    end
end
v_compare =[v_compare;mean(v_compare(1:end-1,:))];
m_compare = [m_compare;mean(m_compare(1:end-1,:))];
p_compare = [p_compare;mean(p_compare(1:end-1,:))];
a_compare = [a_compare;mean(a_compare(1:end-1,:))];

para_est = para_est';
    