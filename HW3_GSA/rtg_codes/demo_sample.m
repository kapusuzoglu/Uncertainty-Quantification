%
% This code samples the scalar QoI and approximates its PDF
%
Ns = 1e4;    
Tf = 80;
y = zeros(Ns, 1);
for i = 1 : Ns
   ratio = 100*(i/Ns); show_progress(ratio);

   xi = 2*rand(3,1)-1;
   y(i) = sir_scalar_qoi(xi, Tf);
end

[f xi] = ksdensity(y);

close all;
figure(1);
hold on;
histnorm(y); 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');


plot(xi, f, '-b', 'linewidth', 2);
xlabel('QoI');
ylabel('distribution');
legend('histogram', 'PDF');
set(gca, 'fontsize', 20);

