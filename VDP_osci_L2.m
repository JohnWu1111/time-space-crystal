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

phi0 = [1;2];
dphi0 = [0;0];

phi = zeros(2,nt);
dphi = zeros(2,nt);
phi(:,1) = phi0;
dphi(:,1) = dphi0;

figure;
x = cospi(phi(:,1));
y = sinpi(phi(:,1));
h = scatter(x,y,'o');
axis([-1 1 -1 1])

for i = 2:nt
    [phi(:,i), dphi(:,i)] = my_runge(phi(:,i-1),dphi(:,i-1),t(i),dt,mu);
    h.XData = cospi(phi(:,i));
    h.YData = sinpi(phi(:,i));
    drawnow
end

phi_diff = phi(1,:) - phi(2,:);

figure;
set(gcf, 'position', [250 70 1500 900]);
titlename = strcat('seed = ',num2str(myseed),'dt = ', num2str(dt));
% subplot(2,1,1)
% plot(t,cos(phi));
% subplot(2,1,2)
% plot(t,cos(phi_diff));
subplot(2,1,1)
plot(t,phi);
subplot(2,1,2)
plot(t,phi_diff);

toc;

function [y, dy] = my_runge(phi, dphi, t, dt, mu)
[c1,d1] = f(phi, dphi, t, mu);
[c2,d2] = f(phi+c1*dt/2, dphi+d1*dt/2, t, mu);
[c3,d3] = f(phi+c2*dt/2, dphi+d2*dt/2, t, mu);
[c4,d4] = f(phi+c3*dt, dphi+d3*dt, t, mu);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
dy = dphi + dt*(d1+2*d2+2*d3+d4)/6;
end

% function [y, dy] = f(x, dx, t, mu)
%     y = dx;
%     dy = -x + mu*(1-x.^2).*dx;
%     diff = x(1)-x(2);
%     dy(1) = dy(1) - 0.5*cos(diff*t);
%     dy(2) = dy(2) + 0.5*cos(diff*t);
% end

function [y, dy] = f(x, dx, t, mu)
    y = dx;
    dy = -x;
    temp = 2*x(1)*x(2);
    dy = dy + cos(temp);
end
