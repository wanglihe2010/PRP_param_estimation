function d = demandcdf(para_est1, para_est2,para_est3, T) 
[time, int, dur, sim_s] = pulse_generation(para_est1(1:4,:),para_est1(5,:),1,1);
demand_test_1 = pulse_aggreation(time, int, dur, sim_s, T);
[f1,x1] = ecdf(demand_test_1(demand_test_1>0));
figure
plot(x1,f1,'r')

[time, int, dur, sim_s] = pulse_generation(para_est2(1:4,:),para_est2(5,:),1,0);
demand_test_2 = pulse_aggreation(time, int, dur, sim_s, T);
[f2,x2] = ecdf(demand_test_2(demand_test_2>0));
hold on
plot(x2,f2,'g')


[time, int, dur, sim_s] = pulse_generation(para_est3(1:4,:),para_est3(5,:),1,0);
demand_test_3 = pulse_aggreation(time, int, dur, sim_s, T);
[f3,x3] = ecdf(demand_test_3(demand_test_3>0));
hold on
plot(x3,f3,'k')
[time, int, dur, sim_s] = pulse_generation(para_est3(1:4,:),para_est3(5,:),2,0);
demand_test_3 = pulse_aggreation(time, int, dur, sim_s, T);
[f3,x3] = ecdf(demand_test_3(demand_test_3>0));
hold on
plot(x3,f3,'k')
[time, int, dur, sim_s] = pulse_generation(para_est3(1:4,:),para_est3(5,:),3,0);
demand_test_3 = pulse_aggreation(time, int, dur, sim_s, T);
[f3,x3] = ecdf(demand_test_3(demand_test_3>0));
hold on
plot(x3,f3,'k')


load('demand_e.mat')

homeint = demand_e{1,3}(:,3);  % unit : L/s
homedur = demand_e{1,3}(:,2);  % unit : s
hometime = demand_e{1,3}(:,1); % unit : s
demand_obs = pulse_aggreation(hometime, homeint, homedur, 214*86400, T);
[f,x] = ecdf(demand_obs(demand_obs>0));
hold on 
plot(x,f)
d=0;