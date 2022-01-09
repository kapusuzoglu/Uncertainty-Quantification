%
% demonstrates impact of fixing parameters
%

Ns = 1e3;
close all

rng(123)
fixing_var_test(2, Ns)
rng(123)
fixing_var_test([2 3], Ns)
rng(123)
fixing_var_test(1, Ns)
