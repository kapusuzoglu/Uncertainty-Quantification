% This code assumes demo_sobol_time has been already run

ndim = 3;
for i = 1 : ndim
   fprintf('param %i analysis \n', i);
   gStot{i} = get_gen_indices(Stot_time{i}, D, timevec);
end

%
% plots
%
close all;
figure(1)
hold on
for i = 1 : ndim
   plot(timevec, gStot{i}, 'linewidth', 2);
end
grid on
legend('p_1', 'p_2', 'p_3');
set(gca, 'fontsize', 20);
xlabel('time');
ylabel('total Sobol indices');
xlim([0 70]);

figure(2);
gStotVec = [gStot{1}(end) gStot{2}(end) gStot{3}(end)];
bar(gStotVec);
set(gca, 'xtick', [1 2 3]);
set(gca, 'xticklabels', {'p_1', 'p_2', 'p_3'});
set(gca, 'fontsize', 20);

