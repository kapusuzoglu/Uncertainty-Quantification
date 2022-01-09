function fixing_var_test(idx, Ns)

for i = 1 : Ns
   xi = 2*rand(3,1)-1;
   xi_fix = xi;
   for k = idx
      xi_fix(k) = 0;
   end
   y(i) = sir_scalar_qoi(xi);
   y_fix(i) = sir_scalar_qoi(xi_fix);
end

[f1 xi1] = ksdensity(y);
[f2 xi2] = ksdensity(y_fix);

figure; hold on;
plot(xi1, f1, 'linewidth',2);
plot(xi2, f2, 'linewidth',2);

param_labels = {'p_1', 'p_2', 'p_3'};

label_str = '';
for k = idx
   label_str = [label_str param_labels{k} ' '];
end
legend('full model', ['model with ' label_str 'fixed']);
xlabel('QoI');
ylabel('PDF');
set(gca, 'fontsize', 20);


