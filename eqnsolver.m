function [param] = eqnsolver(dvar, dmean, dauto, dprob,T)
syms a b t l
    
eqn1 = dmean - T*(a*t*l) == 0;
eqn2 = dvar - (T./t-1+exp(-T./t)) .* (2*l.*(t).^3.*(a.^2+b.^2)) == 0;
eqn3 = dprob - exp(-(t+T)*l)==0;
eqn4 = dauto - 0.5*(1-exp(-1/t*T))^2/(1/t*T-1+exp(-1/t*T)) ==0;
[x] = solve(eqn1,eqn2,eqn3,eqn4);
param =  [ double(x.a); double(x.b); double(x.t); double(x.l)]
