clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

L = 10;
T = 30;
omega = 1;
R = 3;
C = 1;
kappa = 1;
x = (1:L)';
t = 1:T;

pos = x*R - C*cos(kappa*x*R-omega*t);

figure;
h = scatter(t.*ones(L,1),pos,'o');
% axis([-1 1 -1 1])

toc;

