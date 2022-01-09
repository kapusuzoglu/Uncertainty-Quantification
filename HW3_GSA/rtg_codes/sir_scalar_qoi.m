function qoi = sir_scalar_qoi(Xi, Tf)

if nargin < 2
   Tf = 80;
end
% run solver
[sol] = sir_solver(Xi(:), Tf);
timevec = sol.x;
Y = sol.y;
tfinal = sol.x(end);
qoi = (1 / Tf) * trapz(timevec, Y(1,:));
