% UQ
% HW# 1

clc; clear all; close all;


% % Q1
% syms tau w t
% assume(w,'real')
% 
% % R_gg = exp(abs(-tau));
% R_gg = exp(-tau);
% 
% % Eq. 1 from Random Process Simulation tutorial
% S_gg = 1/pi * int(R_gg*cos(w*tau),tau,0,inf);
% 
% 
% N = 500;
% wu = 8*pi;
% dw = wu/N;
% 
% k = 0:1:N-1;
% wk = k*dw;
% 
% w = wk;
% A = sqrt(2*eval(S_gg)*dw);
% 
% % phase angle
% phi = (2*pi) * rand(1,N);
% 
% % gsum = 0;
% % % Eq. 3 from Random Process Simulation tutorial
% % for i=1:N    
% %     gsum = gsum + A(i)*cos(w(i)*t+phi(i));
% % end
% % g = sqrt(2)*gsum;
% 
% 
% % Autocorrelation function of the simulated stochastic process f(t)
% gsum = 0;
% % Eq. 52 from Shinozuka and Deodatis 1991
% for i=1:N    
%     gsum = gsum + A(i)*A(i)*cos(w(i)*t);
% end
% g = gsum;
% 
%  
% 
% % Time step condition
% % dt <= 2*pi/(2*wu) 
% 
% 
% % PLOT
% tau = linspace(0,1,50);
% t = linspace(0,1,50);
% plot(tau,eval(R_gg),t,eval(g),'+')
% xlabel('t');
% ylabel('R(t)');
% legend('Target','Simulation');
% 



% Q5
% syms tau w ws t
% assume(w,'real')
% 
% 
% % R_gg = exp(abs(-tau));
% R_gg = exp(-tau);
% 
% 
% n=9;
% % 
% % 
% % lam = 2/(1+w^2)
% % f = cos(w*t)/sqrt(1+sin(2*w)/2/w)
% % 
% % lam = 2/(1+ws^2)
% % f = sin(ws*t)/sqrt(1-sin(2*ws)/2/ws)
% 
% 
% eq1=w*tan(w)-1==0;
% eq2=w+tan(w)==0;
% sol=solve([eq1,eq2],w)
% sol.w
% % 
% % % zeta
% % zet = randn(1,n);
% % 
% % A = sqrt(lam)*eval(f)*zet;




% Q.5

% 
N=10;
tau = linspace(0,1,N);
rxx =  exp(-abs(-tau));
rr = toeplitz(rxx);
[Psi,D] = eigs(rr,N,'lm');
[lambda,idx] = sort(diag(D),'descend');

Psi = Psi(:,idx);
nn = 10;
n1=100000;
Z = randn(nn,n1);
L = diag( sqrt(lambda(1:nn)) );
X = Psi(:,1:nn)*(L*Z);

% points = linspace(0,1,20)';
R = corr( X' );  % covariance at discrete locations
t = linspace(0,1,N);
plot(tau,rxx)
hold on, plot(t,R(:,1))
xlabel('t');
ylabel('R(t)');
legend('Target','Simulation');
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'KL', 'pdf') %Save figure

% 
% K = @(s,t) exp(-abs(s-t));
% N = 30;
% F = chebop(@(u) fred(K, u));
% [Psi,Lambda] = eigs(F,N,'lm');
% Psi = chebfun(Psi);
% [lambda,idx] = sort(diag(Lambda),'descend');
% Psi = Psi(:,idx);
% 
% % The truncation of the Mercer sum does lead to an underestimate of the 
% % values of the kernel K(s,t). For our example, we should get K(s,s)=1
% % Psi(0,:)*diag(lambda)*Psi(0,:)'
% 
% 
% nn = 10;
% Z = randn(nn,400);
% L = diag( sqrt(lambda(1:nn)) );
% X = Psi(:,1:nn)*(L*Z);
% 
% points = linspace(0,1,N)';
% C = corr( X(points,:)' );  % covariance at discrete locations
% 
% % Plot
% tau = linspace(0,1,N);
% Rgg = exp(-(abs(-tau)));
% plot(tau,Rgg,tau,C(:,1))
% xlabel('t');
% ylabel('R(t)');
% legend('Target','Simulation');
% 
