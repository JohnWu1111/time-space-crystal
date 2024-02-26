clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T = 100;
dt = 1e-3;
t = 0:dt:T;
nt = length(t);
% omega = 1;
mu = 1;

phi0 = [1;1];

phi = zeros(2,nt);
phi(:,1) = phi0;

for i = 2:nt
    phi(:,i) = my_runge(phi(:,i-1),dt,mu);
end

figure;
plot(t,phi(1,:))

toc;

function y = my_runge(phi, dt, mu)
c1 = f(phi, mu);
c2 = f(phi+c1*dt/2, mu);
c3 = f(phi+c2*dt/2, mu);
c4 = f(phi+c3*dt, mu);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = f(x, mu)
    y = x;
    y(1) = x(2);
    y(2) = -x(1) + mu*(1-x(1)^2)*x(2);
end
