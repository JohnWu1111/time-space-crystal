clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

T = 10;
dt = 1e-3;
t = 0:dt:T;
nt = length(t);
omega = 1;
J = 2;

phi0 = 2*rand(3,1);

phi = zeros(3,nt);
phi(:,1) = phi0;

phi_diff = zeros(3,nt);
phi_diff(:,1) = circshift(phi(:,1),1) - phi(:,1);

x = cospi(phi(:,1));
y = sinpi(phi(:,1));
% figure;
% h = scatter(x,y,'o');
% axis([-1 1 -1 1])

for i = 2:nt
    dphi = my_runge(phi(:,i-1),dt,omega,J);
    phi(:,i) = phi(:,i-1) + dt*dphi;
    phi_diff(:,i) = circshift(phi(:,i),1) - phi(:,i);
%     h.XData = cospi(phi(:,i));
%     h.YData = sinpi(phi(:,i));
%     drawnow
end

order = cospi(phi);
order_ST = sum(abs(order));
phi_diff = mod(phi_diff+1,2)-1;

figure;
set(gcf, 'position', [250 70 1500 900]);
titlename = strcat('seed = ',num2str(myseed),'dt = ', num2str(dt), 'omega = ', num2str(omega));
subplot(3,1,1)
plot(t,order');
subplot(3,1,2)
plot(t,phi_diff');
subplot(3,1,3)
plot(t,order_ST)

toc;

function y = my_runge(phi, dt, omega, J)
c1 = f2(phi, omega, J);
c2 = f2(phi+c1*dt/2, omega, J);
c3 = f2(phi+c2*dt/2, omega, J);
c4 = f2(phi+c3*dt, omega, J);
y = (c1+2*c2+2*c3+c4)/6;
end

function y = f2(x, omega, J)
%     y = x;
%     y(1) = omega + cos((x(1)-x(2))/2);
%     y(2) = omega - cos((x(1)-x(2))/2);
    x_diff = mod(x-x'+1,2)-1;
    for i = 1:3
        x_diff(i,i) = 1;
    end
    y = omega + J*sum(1./x_diff,2) -1;
end

function y = f1(x, omega, J)
%     y = x;
%     y(1) = omega + cos((x(1)-x(2))/2);
%     y(2) = omega - cos((x(1)-x(2))/2);
    x_diff = mod(x-x'+1,2)-1;
    fact = sign(x_diff);
    diff = pot(x_diff);
    y = omega + J*sum(fact.*diff,2);
end

function y = pot(x)
y = x;
for i = 1:3
    for j = 1:3
        if abs(x(i,j)) < 2/3
            y(i,j) = cospi(3*x(i,j)/4);
        else
            y(i,j) = 0;
        end
    end
end
end
