%
% solves the system at a nominal parameter and shows solution
%

% nominal parameter:
x0 = [0 0 0]';

Tf = 80;
[sol] = sir_solver(x0, Tf);

%
% postprocess
%

% these are the nominal parameters  
params = get_physical_params(x0);
p1 = params(1);   
p2 = params(2);
p3 = params(3);
 
p1_str = num2str(p1, '%g');
p2_str = num2str(p2, '%g');
p3_str = num2str(p3, '%g');

figure(1) 
plot(sol.x, sol.y, 'linewidth', 2);
legend('S', 'I', 'R');
title(['p_1 = ' p1_str ', p_2 = ' p2_str ', p_3 = ' p3_str]); 
xlabel('time');
ylabel('population');
set(gca, 'fontsize', 20);



