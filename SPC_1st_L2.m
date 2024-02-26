clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T = 30;
dt = 1e-3;
t = 0:dt:T;
nt = length(t);
omega1 = 0;
omega2 = 1;

phi0 = 2*rand(2,1);

phi = zeros(2,nt);
phi(:,1) = phi0;

x = cospi(phi(:,1));
y = sinpi(phi(:,1));
% figure;
% h = scatter(x,y,'o');
% axis([-1 1 -1 1])

for i = 2:nt
    phi(:,i) = my_runge(phi(:,i-1),t(i-1),dt,omega1,omega2);
%     h.XData = cospi(phi(:,i));
%     h.YData = sinpi(phi(:,i));
%     drawnow
end

phi_diff = phi(1,:) - phi(2,:);

figure;
set(gcf, 'position', [250 70 1500 900]);
titlename = strcat('seed = ',num2str(myseed),'dt = ', num2str(dt));
subplot(3,1,1)
plot(t,phi);
subplot(3,1,2)
plot(t,mod(phi-omega1*t,2),t,mod(phi(2,:)-omega1*t+1,2));
% plot(t,phi,t,phi(2,:)-1);
subplot(3,1,3)
plot(t,phi_diff);

toc;

function y = my_runge(phi, t, dt, omega1,omega2)
c1 = f(phi, t, omega1,omega2);
c2 = f(phi+c1*dt/2, t+dt/2, omega1,omega2);
c3 = f(phi+c2*dt/2, t+dt/2, omega1,omega2);
c4 = f(phi+c3*dt, t+dt, omega1,omega2);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = f(x, t, omega1,omega2)
    y = x;
    diff = x(1)-x(2);
    y(1) = omega1 + cospi(omega2*t) - cospi(diff/2);
    y(2) = omega1 - cospi(omega2*t) + cospi(diff/2);
end
