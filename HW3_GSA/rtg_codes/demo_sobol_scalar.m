%
% Computes the Sobol' indices for the scalar QoI
%
ndim = 3;
Ns = 1e3;

F = @sir_scalar_qoi;
for i = 1 : ndim
   fprintf('param %i analysis \n', i);
   [Xi F1, F2, mu, D, Stot(i)] = get_sobol_scalar(F, ndim, Ns, i);
end

close all;
figure(1);
param_labels = {'p_1', 'p_2', 'p_3'}; 
bar(Stot);
set(gca,'fontsize', 20, 'xticklabels', param_labels); 
xlim([0 4])
ylim([0 1])
ylabel('total Sobol index')
