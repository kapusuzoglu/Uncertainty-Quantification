clc;
clear;
close all
rng(1);

N = 16;
n = N^3-5*N;


%% Assume uncorrelated inputs
% Grid Samples
% P~Type I
% Type I Extreme Value Distribution: 
% CDF: F_Y = exp(-exp(-alphan*(y-u)))
% PDF: f_Y = alphan*exp(-alphan*(y-u))*exp(-exp(-alphan*(y-u)))
% alphan = 1/sqrt(6)*(pi/sigmaYn)
% u = mu_Y - 0.5772/alphan
mu_x1 = 2000;
cov_x1 = 0.2;
sigma_x1 = mu_x1*cov_x1;

alphan = 1/sqrt(6)*(pi/sigma_x1);
u = mu_x1 - 0.5772/alphan;
F = exp(-exp(-alphan*(mu_x1-u)));
f = alphan*exp(-alphan*(mu_x1-u))*exp(-exp(-alphan*(mu_x1-u)));
sigmaN_x1 = normpdf(norminv(F))/f;
muN_x1 = mu_x1 - sigmaN_x1*(norminv(F));

% E~Lognormal
mu_x2 = 30000;
cov_x2 = 0.1;
sigma_x2 = mu_x2*cov_x2;
zeta_x2 = sqrt(log(1+cov_x2^2));
lambda_x2 = log(mu_x2) - 0.5*zeta_x2^2;
% sigmaN_x2 = zeta_x2*mu_x2
% muN_x2 = mu_x2*(1-log(mu_x2)+lambda_x2)
muN_x2 = log(mu_x2 / sqrt(1 + sigma_x2^2/mu_x2^2));
sigmaN_x2 = log(1 + sigma_x2^2/mu_x2^2);

% I~Normal
muN_x3 = 10;
cov_x3 = 0.05;
sigmaN_x3 = muN_x3*cov_x3;

dev = 3;
P = linspace( mu_x1 - dev*sigma_x1, mu_x1 + dev*sigma_x1, N);
E = linspace( mu_x2 - dev*sigma_x2, mu_x2 + dev*sigma_x2, N);
I = linspace( muN_x3 - dev*sigmaN_x3, muN_x3 + dev*sigmaN_x3, N);
[PP, EE, II] = ndgrid(P,E,I);
% Input

Nsamples = 4096;
% Uniform numbers
r = rand(Nsamples,3);
x1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1);
x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2);
x3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3);

input_sample1 = [x1 x2 x3];

Y = func_beam(input_sample1);

tic
%% Double loop method
First_order_double_loop   = DoubleLoopGSA_First(input_sample1, @func_beam, n);
Total_effects_double_loop = DoubleLoopGSA_Total(input_sample1, @func_beam, n);
toc


% dev = 3;
% P = linspace( mu_x1 - dev*sigma_x1, mu_x1 + dev*sigma_x1, N);
% E = linspace( mu_x2 - dev*sigma_x2, mu_x2 + dev*sigma_x2, N);
% I = linspace( muN_x3 - dev*sigmaN_x3, muN_x3 + dev*sigmaN_x3, N);
% [PP, EE, II] = ndgrid(P,E,I);
% % Input
% input_sample1 = [PP(:) EE(:) II(:)];
% 
% Y = func_beam(input_sample1);

tic
%% Single loop methods
% CDF
input_cdf = {[1 2 3], {'Extreme Value','Lognormal','Normal'}, ...
    {mu_x1,sigma_x1}, {lambda_x2,zeta_x2}, {muN_x3,sigmaN_x3}};
SingleFirst_Index = SingleLoop_FirstOrder(input_sample1, ...
    Y, sqrt(N)+N,  input_cdf);
toc


SingleFirst_Index_Algo1 = SingleLoop_FirstOrder_Algo1(input_sample1, ...
    Y, sqrt(N)+N,  input_cdf);



N = 2e3;
n = ceil(sqrt(N));
%% Latin Hypercube Samples
Mu = [muN_x1,muN_x2,muN_x3]; % mean value
Std = Mu.*[cov_x1 cov_x2 cov_x3]; % standard deviation
Sigma = diag(Std.^2); % standard deviation;
input_sample2 = lhsnorm(Mu, Sigma, N);
A = input_sample2;
B = lhsnorm(Mu, Sigma, N);

tic
First_order_Saltelli = Sen_FirstOrder_Saltelli(A, B, @func_beam);
Total_effects_Saltelli = Sen_TotalEffect_Saltelli(A, B, @func_beam);
toc

%% PCE

% Input in terms of Standard Random Variates
x1 = icdf('Normal',cdf('Extreme Value', input_sample1(:,1), mu_x1, sigma_x1), 0, 1);
x2 = (log(input_sample1(:,2))-lambda_x2)/zeta_x2;
x3 = ((input_sample1(:,3)-muN_x3)/sigmaN_x3);

%% Linear regression analysis to find the trend
X = [ones(size(Y)) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3];
b = regress(Y, X);

%% MCS on the original model
%--------------------------------------------------------------------------
N = 1000;
% Uniform numbers
r = rand(N,3);
x1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1);
x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2);
x3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3);
% MCS_disp = x1*L./(x2.*x3);
%-----------------------------------------------------------------------------------------------------
% End of MCS

% Surrogate model - PCE
x1 = icdf('Normal',cdf('Extreme Value', x1, mu_x1, sigma_x1), 0, 1);
x2 = (log(x2)-lambda_x2)/zeta_x2;
x3 = ((x3-muN_x3)/sigmaN_x3);

Y_pce = sum(b'.*[ones(N,1) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3],2);

psi = [ones(N,1) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3];
psi2 = psi.^2;
b2 = b.^2;

% VY = var(Y_pce)

VY = b2(2:end)'*mean(psi2(:,2:end))';

% c*psi_1
cpsi1 = [ b2(2) b2(5)]*mean([psi2(:,2) psi2(:,5)])';
cpsi2 = [ b2(3) b2(6)]*mean([psi2(:,3) psi2(:,6)])';
cpsi3 = [ b2(4) b2(7)]*mean([psi2(:,4) psi2(:,7)])';

Si_1 = cpsi1/VY
Si_2 = cpsi2/VY
Si_3 = cpsi3/VY
Si_1+Si_2+Si_3
FirstIndex_PCE=[Si_1 Si_2 Si_3];

cpst1 = [ b2(2) b2(5) b2(8) b2(9)]*mean([psi2(:,2) psi2(:,5) psi2(:,8) psi2(:,9)])';
cpst2 = [ b2(3) b2(6) b2(8) b2(10)]*mean([psi2(:,3) psi2(:,6) psi2(:,8) psi2(:,10)])';
cpst3 = [ b2(4) b2(7) b2(9) b2(10)]*mean([psi2(:,4) psi2(:,7) psi2(:,9) psi2(:,10)])';

ST_1 = cpst1/VY;
ST_2 = cpst2/VY;
ST_3 = cpst3/VY;
TotalIndex_PCE=[ST_1 ST_2 ST_3];
ST_1+ST_2+ST_3







% % Sobol' method
% First_order_Sobol = Sen_FirstOrder_Sobol_07(A, B, @func_beam);
% 
% % improved FAST based on RBD
% k = size(Mu, 2);
% First_order_FAST = Sen_FirstOrder_RBD(N, k, @beam_RBD, Mu, Std);

%% plot
figure(1)
plot(First_order_double_loop, 'r+'); hold on;
plot(SingleFirst_Index_Algo1, 'm*'); hold on;
plot(SingleFirst_Index, 'bs'); hold on;
plot(FirstIndex_PCE, 'go'); hold on;
plot(First_order_Saltelli, 'd'); hold on;
% plot(First_order_Sobol, 'k^'); hold on;
% plot(First_order_FAST, 'mv');
xlim([0, 4]);
ylim([0, 1]);
set(gca, 'XTickLabel',{'P','E','I'},'XTick',[1 2 3]);
xlabel('input')
ylabel('First-order index')
% legend({'Double loop', 'Saltelli', 'Sobol', 'FAST'}, 'location', 'northwest')
legend({'Double loop Original Model','Single loop Algorithm 1', 'Single loop Algorithm 2', 'PCE','Latin Hypercube Sampling'}, 'location', 'northeast')
title('First-order index, uncorrelated input')

figure(2)
plot(Total_effects_double_loop, 'r+'); hold on;
plot(TotalIndex_PCE, 'go'); hold on;
plot(Total_effects_Saltelli, 'bs');
xlim([0, 4]);
ylim([0, 1]);
set(gca, 'XTickLabel',{'P','E','I'},'XTick',[1 2 3]);
xlabel('input')
ylabel('Total effects index')
legend({'Double loop', 'PCE','Latin Hypercube Sampling'}, 'location', 'northeast')
title('Total effects index, uncorrelated input')

% %% Assume correlated inputs
% corr = [1	0.174	0.451	0.082	-0.134	0.004;
% 0.174	1	-0.8	0.059	-0.125	-0.082;
% 0.451	-0.8	1	-0.004	0.033	0.08;
% 0.082	0.059	-0.004	1	-0.105	-0.4;
% -0.134	-0.125	0.033	-0.105	1	0.279;
% 0.004	-0.082	0.08	-0.4	0.279	1]; % correlation matrix
% 
% COV = zeros(k, k);
% for i = 1:k
%     for j = 1:k
%         COV(i, j) = corr(i, j) * Std(i) * Std(j);
%     end
% end
% First_order_double_loop_2 = GSA_FirstOrder_mvn(Mu', COV, @func_beam, 1000);
% 
% figure(3)
% plot(First_order_double_loop_2, 'bs'); hold on;
% xlim([0, 7]);
% set(gca, 'XTickLabel',{'P','E','v','b','H','L'},'XTick',[1 2 3 4 5 6]);
% xlabel('input')
% ylabel('First-order index')
% legend({'Double loop'}, 'location', 'northwest')
% title('Total effects index, correlated input')