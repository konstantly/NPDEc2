function [T,Q,P,H] = HamSolver(q0, p0, N, h, Force, Method, Pars)
% HamSolver solves the Hamiltonian System IVP 
%    dq/dt = p; dp/dt = F(q), 
% q(0)=q0, p(0)=p0.  
%
% The force F(q)=-nabla U(q), for potential energy U(q).
%
% The method takes N steps of size h with an ODE 
% numerical method idenfified by the Method argument
%
% The Pars structure contains parameters which are passed to the 
% routine Force.  
% Pars.Mass is a column vector of masses, so that M = diag(Pars.Mass) is the mass matrix
% In the plane the masses should be repeated once, i.e. 
% Pars.Mass = [ m1 m1 m2 m2...]
%
% Force identifies a function that takes two input arguments: q Pars.
% and has 2 output arguments:  -nabla U and U (potential energy).
%
% HamSolver Returns a vector of timesteps T and the arrays Q,P whose columns represent
% the position vectors and momentum vectors at successive 
% time points T=[0,h, 2h,...Nh].  It also returns H, a vector of energies H(q,p)=p^Tp/2 + U(q)
% at the successive time points.
%
n = length(q0); % q0 is assumed to be a column vector, n is 
                % number of degrees of freedom
%
Q=zeros(n,N+1); % array of position vectors
P=zeros(n,N+1); % array of momentum vectors
T=zeros(1,N+1); % array of time values
H=zeros(1,N+1); % array of energies

global F;
global U; 
% we use global variables to keep track of the most recent
% evaluations of the force and potential ; allows to avoid 
% duplication of force evaluations.

% initial values
Q(:,1) = q0;
P(:,1) = p0;
T(1) = 0;
[F,U] = Force(q0, Pars); % the force function should return both force and potential
H(1) = sum(p0.*(p0./Pars.Mass))/2 + U;

% loop over timesteps, computing successive points along solution
for n=1:N
  [Q(:,n+1),P(:,n+1)] = ODEStep(Q(:,n),P(:,n),h,Force,  Method, Pars);
  T(n+1) = T(n)+h;
  H(n+1) = sum(P(:,n+1).*(P(:,n+1)./Pars.Mass))/2 + U;
end
%
%
function [qnp1, pnp1] = ODEStep(qn,pn, h,Force,  Method, Pars)
%
% Takes a single step using a defined ODE method for solving the Hamiltonian
% system (potential energy and force defined by Force) with stepsize h. 
%
% Options for Method include 'Euler', 'RK4' and 'Verlet'
% 
% The force field defined by Force should be the name of an m-file in single quotes
% or else anonymous function name. Force should take the arguments q and Pars
% where Pars is a parameter structure whose elements can be any matlab
% object. Pars.Mass is a column vector of masses.  After qnp1 is computed, the force should also be
% computed and stored in the global variable F (and the energy in U)
% to get ready for the next step and for energy recording.
% 
%
global F;
global U;

switch(Method)
    case {'Euler', 'euler'}
        qnp1 = qn + h*(pn./Pars.Mass);
        pnp1 = pn + h*F;
        [F,U] = Force(qn, Pars); %evaulate the force and potential
    case {'RK4', 'rk4'}
        P1 = pn;
        Q1 = qn;
        F1 = F;
        P2 = pn + (h/2)*F1;
        Q2 = qn + (h/2)*(P1./Pars.Mass);
        F2 = Force(Q2, Pars);
        P3 = pn + (h/2)*F2;
        Q3 = qn + (h/2)*(P2./Pars.Mass);
        F3 = Force(Q3, Pars);
        P4 = pn + (h)*F3;
        Q4 = qn + (h)*(P3./Pars.Mass);
        F4 = Force(Q4, Pars);
        qnp1 = qn + (h/6)*((P1 + 2*P2 + 2*P3 + P4)./Pars.Mass);
        pnp1 = pn + (h/6)*(F1 + 2*F2+2*F3+F4);
        [F,U] = Force(qnp1, Pars);
    case {'Verlet', 'leapfrog'}
        pnph = pn + .5*h*F;
        qnp1 = qn + h*(pnph./Pars.Mass);
        [F,U] = Force(qnp1, Pars);
        pnp1 = pnph + .5*h*F;
    otherwise,
        disp('Method not implemented')
end
