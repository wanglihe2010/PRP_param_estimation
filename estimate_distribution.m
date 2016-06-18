load('inputparam.mat');
n_mc = 500;
tau = zeros(n_mc,24);
tau_n = zeros(n_mc,24);
for i = 1:n_mc

[time, int, dur, sim_s] = pulse_generation(inputparam(1:3,:),inputparam(4,:), 1, 0);
times = 2;
t1 = 60;
demand1 = pulse_aggreation(time, int, dur, sim_s, t1);
demand2 = sum(reshape(demand1,times, length(demand1)/times)',2);
[para_est] = para_estimation_24(demand1,demand2,t1,times*t1)';

[para_est_n] = para_estimation_new(demand1,t1);

tau(i,:) = para_est(3,:);
tau_n(i,:) = para_est_n(3,:);

i
end

tau_error = abs(tau - ones(100,1)*inputparam(3,:))./(ones(100,1)*inputparam(3,:))*100;
tau_n_error = abs(tau_n - ones(100,1)*inputparam(3,:))./(ones(100,1)*inputparam(3,:))*100;