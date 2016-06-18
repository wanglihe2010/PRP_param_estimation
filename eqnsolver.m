function [tau] = eqnsolver( dauto,T)
syms t 
    
eqn4 = dauto - 0.5*(1-exp(-1/t*T))^2/(1/t*T-1+exp(-1/t*T)) ==0;
[x] = solve(eqn4);
tau =  [double(x.t)];
