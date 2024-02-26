clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

T = 50;
dt = 1e-3;
t = 0:dt:T;
nt = length(t);
omega1 = 2;
omega2 = 1;

phi0 = 2*pi*rand(2,1);
dphi0 = 2*pi*rand(2,1);

phi = zeros(2,nt);
phi(:,1) = phi0;
dphi = zeros(2,nt);
dphi(:,1) = dphi0;

x = cos(phi(:,1));
y = sin(phi(:,1));
% figure;
% h = scatter(x,y,'o');
% axis([-1 1 -1 1])

for i = 2:nt
    [phi(:,i),dphi(:,i)] = my_runge(phi(:,i-1),dphi(:,i-1),t(i-1),dt,omega1,omega2);
%     h.XData = cos(phi(:,i));
%     h.YData = sin(phi(:,i));
%     drawnow
end

phi_diff = phi(1,:) - phi(2,:);

figure;
set(gcf, 'position', [250 70 1500 900]);
titlename = strcat('seed = ',num2str(myseed),'dt = ', num2str(dt));
subplot(3,1,1)
plot(t,phi);
subplot(3,1,2)
% plot(t,mod(phi-omega1*t,2),t,mod(phi(2,:)-omega1*t+1,2));
plot(t,phi(1,:),t,circshift(phi(2,:),1/dt));
subplot(3,1,3)
plot(t,phi_diff);

toc;

function [y, dy] = my_runge(phi, dphi, t, dt, omega1,omega2)
[c1,d1] = f(phi, dphi, t, omega1,omega2);
[c2,d2] = f(phi+c1*dt/2, dphi+d1*dt/2, t, omega1,omega2);
[c3,d3] = f(phi+c2*dt/2, dphi+d2*dt/2, t, omega1,omega2);
[c4,d4] = f(phi+c3*dt, dphi+d3*dt, t, omega1,omega2);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
dy = dphi + dt*(d1+2*d2+2*d3+d4)/6;
end

function [y, dy] = f(x, dx, t, omega1,omega2)
    y = dx;
    dy = x;
    diff = x(1)-x(2);
%     dy(1) =  - omega1*diff;
%     dy(2) =  - omega2*diff;
    dy(1) = (1-x(1)^2)*x(1) - 2*omega1*cos(diff/2);
    dy(2) = (1-x(2)^2)*x(2) + 2*omega1*cos(diff/2);
end
