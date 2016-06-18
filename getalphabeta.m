function [alpha_beta] = getalphabeta(dmean, dvar, tau, lambda, T)
option = 1; % 1 - equation solution
            % 2 - least sqaure solution
            % 3 - solve func solution
if option == 1
    alpha = dmean / tau/lambda/T;
    beta = sqrt(dvar/(2*lambda*(tau)^3*(T/tau-1+exp(-T/tau)))- alpha^2); 
    if ~isreal(beta)
        option = 2;
    end
end
if option == 2
    l = leastsquare([0.2, 0.1], dvar,dvar,dmean,dmean,tau,lambda,T,T);
    if ~isfinite(l)
        c = 0;
    end
    options = optimset('Algorithm','interior-point','Hessian','bfgs','Display','off');
    [x,fval,existflag] = fmincon('leastsquare',[0.2,0.1],[],[],[],[],[0 0 ],[],[],options,dvar,dvar,dmean,dmean,tau,lambda,T,T);

    alpha = x(1);
    beta = x(2);
    l = leastsquare(x, dvar,dvar,dmean,dmean,tau,lambda,T,T)
end
if option == 3
    syms a b 
    
    eqn1 = dmean - T*(a*tau*lambda) == 0;
    eqn2 = dvar - (T./tau-1+exp(-T./tau)) .* (2*lambda.*(tau).^3.*(a.^2+b.^2)) == 0;
    [x] = solve(eqn1,eqn2);
    [alpha, beta] =  [ double(x.a) double(x.b)];
end
alpha_beta = [alpha beta];