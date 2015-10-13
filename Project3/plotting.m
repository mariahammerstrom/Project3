% Plotting of integrand

n = 100;                      % no. of points
a = 0;                      % lower limit
b = 10;                     % upper limit
alpha = 2.0;

x = linspace(a,b,n);
f = exp(-alpha*x);        % function
exp(-2*3)

plot(x,f)
title('Visualization of integrand')
xlabel('$r_i$','Interpreter','latex')
ylabel('$f(r_i)$','Interpreter','latex')