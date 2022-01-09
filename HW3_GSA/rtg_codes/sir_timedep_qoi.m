function qoi = sir_timedep_qoi(Xi, timevec)

% run solver
[sol] = sir_solver(Xi, timevec(end));
spec_idx = 2;
qoi = deval(sol,timevec, spec_idx);

