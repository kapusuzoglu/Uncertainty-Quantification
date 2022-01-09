function [sol] = sir_solver(xi, tfinal)
%
%   xi is a 3 x 1 vector with entries in [-1, 1]
%

S0 = 499;
I0 = 1;
R0 = 0;
Y0 = [S0; I0; R0];


% map xi entries to physical ranges
params = get_physical_params(xi);

if nargin < 2
   tfinal = 100;
end
tol = 1e-7;
ode_options = odeset('RelTol', tol);

frhs = @(t,y)(SIR_rhs(t, y, params));
sol = ode45(frhs, [0 tfinal], Y0, ode_options);

% rhs of ODE system
function dy = SIR_rhs(t,y,params);
p1 = params(1);   
p2 = params(2);
p3 = params(3);

beta = p1 + p2 * sin(t);
dy = [-beta*y(2)*y(1);
      beta*y(2)*y(1) - p3 * y(2);
      p3*y(2) 
     ];
