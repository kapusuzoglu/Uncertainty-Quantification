
%
% Computes the DGSMs and computes bounds on Sobol' indices
% (assumes demo_sobol_scalar has been already been run)
%

Ns = 500;
ep = 1e-4;
ndim = 3;
E = eye(ndim);
nu = 0;
y = zeros(Ns, 1);
dy = zeros(Ns, 3);

f = @(x)(sir_scalar_qoi(x));
for i = 1 : Ns
   ratio = 100*(i/Ns);
   show_progress(ratio);

   xi = 2*rand(3,1)-1;
   % function evaluation
   y(i) = f(xi);
   
   % gradient (with FD)
   for k = 1 : ndim
       yp(k) = f(xi + ep*E(:, k));
       dy(i, k) = (yp(k) - y(i)) / ep;
   end
   nu = nu + dy(i,:).^2;
end
nu = nu / Ns;
D = var(y(:));
bound = 4/pi^2 * nu/D;

figure(2);
bar([Stot(:) bound(:)]);
set(gca, 'fontsize', 20);
legend('S_i^{tot}', 'Bound');
param_labels = {'p_1', 'p_2', 'p_3'};
set(gca,'fontsize', 20, 'xticklabels', param_labels);

%figure(3);
%subplot(3,1,1)
%histnorm(dy(:,1).^2);
%subplot(3,1,2)
%histnorm(dy(:,2).^2);
%subplot(3,1,3)
%histnorm(dy(:,3).^2);

