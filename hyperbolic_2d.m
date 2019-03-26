clear all;
close all;

% Define a "spatial" grid
M = 100;
D = 1;
x = linspace(0, D, M + 1);
y = linspace(0, D, M + 1);
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';
% X(m, p) = x_{m - 1} for all p
% Y(m, p) = y_{p - 1} for all m
dx = D / M;

chi = chi(x,y);

% Define information about the step sizes in the "time" dimension
N = 100;
T = 100;
dt = 0.1;

% The initial condition
psi_nm1 = Pxy(x,y);



% First "timestep"
psi_h = 
psi_n = 

% Remaining "timesteps"
for n = 1:N - 1
  psi_np1 = psi_nm1 - 2*dt*J
  
  psi_nm1 = psi_n;
  psi_n = psi_np1;
end

% For two-dimensional plotting you may wish to use contourf or pcolor
