function [Xi F1, F2, mu sigma2 Sy_tot Sy] = get_sobol_scala(F, ndim, Ns, yidx)
%
%  F    : function handle to simulation routine
%         the goal is to assess sensitivity of F(x_1, ..., x_ndim) to parameters, x_1, ..., x_ndim
%
%  ndim : stochastic dimension 
%
%  Ns   : the number of samples used to estimate the sensitivity indices
%
%  yidx : the index set for which we want to compute the index,
%
%         example yidx = [1 2] gives S_12, yidx = [2] gives S_2, etc.
%         the following is formulated following the paper of Sobol ...
% 
% reference: 
%   Sobol, I.M., Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates, 
%   Mathematics and Computers in Simulation, 2001.

m = length(yidx);
p = ndim - m;

% get the basic samples needed
eta        = 2*rand(Ns, m)-1;
etaprime   = 2*rand(Ns, m)-1;
zeta       = 2*rand(Ns, p)-1;
zetaprime  = 2*rand(Ns, p)-1;

% index sets
I = [1 : ndim];     
idx = setdiff(I, yidx);

% set up the samples
Xi1 = zeros(Ns, ndim);
Xi2 = zeros(Ns, ndim);
Xi3 = zeros(Ns, ndim);

Xi1(:, yidx) = eta;
Xi1(:, idx)  = zeta;

Xi2(:, yidx)   = eta;
Xi2(:, idx) = zetaprime;

Xi3(:, yidx)   = etaprime;
Xi3(:, idx) = zeta;

% get the needed function evaluations 
F1 = zeros(Ns, 1);
F2 = zeros(Ns, 1);
F3 = zeros(Ns, 1);
% the sampling loop
for k = 1 : Ns
   ratio = (k/Ns)*100;
   show_progress(ratio);

   F1(k,:) = F(Xi1(k, :));
   F2(k,:) = F(Xi2(k, :));
   F3(k,:) = F(Xi3(k, :));
end

% Compute the estimators
phi  = F1;
phi2 = F1.^2;
psi  = phi .* F2;
chi  = 0.5 * (phi - F3).^2;

% compute the results
mu      = mean(phi);
mu2     = mean(phi2);
sigma2  = mu2 - mu.^2;
mupsi   = mean(psi);
Sy      = (mupsi - mu.^2) ./ sigma2;
Sy_tot  = mean(chi) ./ sigma2;

Xi = Xi1;

mu     = mu(:);
sigma2 = sigma2(:);
Sy     = Sy(:);
Sy_tot = Sy_tot(:);
