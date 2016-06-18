function [meanobs, varobs, Pobs, autoobs] = statistics_timescale_obs(hometime, homeint, homedur, totaltime_s, aggr_time)
% observed second order statitics on different aggregation time 
% n*24
if nargin <5
    
    for i = 1:3600
        if mod(3600,i) ==0
            b(i) = 1;
        end
    end
    aggr_time  = find(b == 1)';
end


n = length(aggr_time);
Pobs = zeros(n,24);
varobs = zeros(n,24);
meanobs = zeros(n,24);
autoobs = zeros(n,24);
for i = 1:n
    t1 =aggr_time(i);
    homedemand = pulse_aggreation(hometime, homeint, homedur, totaltime_s, t1);
    demand_reshape1 = reshape(homedemand, 24*3600/t1, length(homedemand)/24/3600*t1);
    demand_reshape= zeros(214*3600/t1,24);
    for ii = 1:24
    demand_reshape(:,ii) = reshape(demand_reshape1(3600/t1*(ii-1)+1:3600/t1*ii,:),214*3600/t1,1);
    end
    
    %demand_reshape = reshape(demand_reshape1',length(homedemand)/24,24);

   
     for j = 1:24
        Pobs(i,j) = length(find(demand_reshape(:,j)==0))/length(demand_reshape(:,j));
        a = autocorr(demand_reshape(:,j),1);
        autoobs(i,j) = a(2);
     end
     meanobs(i,:) = mean(demand_reshape);
     varobs(i,:) = var(demand_reshape);
end
