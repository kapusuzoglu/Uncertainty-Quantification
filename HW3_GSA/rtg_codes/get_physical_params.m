function params = get_physical_params(x)

% map x entries to physical ranges
meanpar = [0.001; .002; .1];
a = meanpar - 0.2 * meanpar;
b = meanpar + 0.2 * meanpar;
params  = 0.5 * (a + b) + 0.5 * (b - a) .* x;
