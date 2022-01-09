timevec = [0.01 : 0.1 : 70];
ndim = 3;
Ns = 1e3;
Nt = length(timevec);

F = @(Xi)(sir_timedep_qoi(Xi(:), timevec));
for i = 1 : ndim
   fprintf('param %i analysis \n', i);
   [Xi F1, F2, mu, D, Stot_time{i}] = get_sobol_time(F, ndim, Nt, Ns, i);
end

%
% plots
%
close all;
figure(1)
hold on
for i = 1 : ndim
   plot(timevec, Stot_time{i}, 'linewidth', 2);
end
grid on
legend('p_1', 'p_2', 'p_3');
set(gca, 'fontsize', 20);
xlabel('time');
ylabel('total Sobol indices');
xlim([0 70]);

figure(2);
hold on;
plot(timevec, mu, '-b', 'linewidth', 2);
plot(timevec, mu-2*sqrt(D), 'r--', 'linewidth', 2);
plot(timevec, mu+2*sqrt(D), 'r--', 'linewidth', 2);
legend('mean', '\pm 2 std dev');
xlabel('time');
ylabel('infected population');
set(gca, 'fontsize', 20);
xlim([0 70]); 

figure(3);
hold on;
plot(timevec, sqrt(D), 'linewidth', 2);
xlabel('time');
ylabel('std deviation');
set(gca, 'fontsize', 20);
xlim([0 70]);
