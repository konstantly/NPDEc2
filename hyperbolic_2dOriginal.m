clear all;
close all;

% Define a "spatial" grid
M = 
D = 1;
x = linspace(0, D, M + 1);
y = linspace(0, D, M + 1);
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';
% X(m, p) = x_{m - 1} for all p
% Y(m, p) = y_{p - 1} for all m
dx = D / M;

chi = 

% Define information about the step sizes in the "time" dimension
N = 
T = 
dt = 

% The initial condition
psi_nm1 = 

% First "timestep"
psi_h =
psi_n =

% Remaining "timesteps"
for n = 1:N - 1
  psi_np1 =
  
  psi_nm1 = psi_n;
  psi_n = psi_np1;
end

% For two-dimensional plotting you may wish to use contourf or pcolor
