% Plotting of integrand

n = 100;                      % no. of points
a = 0;                      % lower limit
b = 10;                     % upper limit
alpha = 2.0;

x = linspace(a,b,n);
f = exp(-2*alpha*x);        % function

plot(x,f)
title('Visualization of integrand')
xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$','Interpreter','latex')