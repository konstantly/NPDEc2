%% a
%[a, b, c] = Pxy(1,2,@Jxp)

x = [1, 2, 3];
y = [1, 2, 2];

[X, Y] = meshgrid(x, y);
X = X';
Y = Y';
yinit = sin(pi*X).*sin(pi*Y)

chi = sin(pi*X).*sin(pi*Y)



%% b

clear all;
close all;
%%

% Define a "spatial" grid
M = 10;
D = 1;
x = linspace(0, D, M + 1);
y = linspace(0, D, M + 1);
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';
% X(m, p) = x_{m - 1} for all p
% Y(m, p) = y_{p - 1} for all m
dx = D / M;

chi = sin(pi*X).*sin(pi*Y);

% Define information about the step sizes in the "time" dimension
N = 100;
T = 100;
dt = 0.1;

% The initial condition
psi_nm1 = exp(-((x-0.25).^2+(y-0.6).^2)/0.08).*chi;









