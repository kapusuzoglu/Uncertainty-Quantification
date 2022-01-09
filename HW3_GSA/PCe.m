% THE CANTILEVER BEAM SUBJECTED TO CONCENTRATED LOAD
%% METHODS: POLYNOMIAL CHAOS, MONTE CARLO SIMULATION & GAUSSIAN PROCESS MODEL
% x1=P, x2=E and x3=I are random fields
clc;
clear all;close all;
%--------------------------------------------------------------------------
%% Polynomial Chaos Expansion Model
%--------------------------------------------------------------------------
%collocation points
P = [-1.73205080756888,-1.73205080756888,-1.73205080756888;
    0,-1.73205080756888,-1.73205080756888;
    1.73205080756888,-1.73205080756888,-1.73205080756888;
    -1.73205080756888,0,-1.73205080756888;
    0,0,-1.73205080756888;
    1.73205080756888,0,-1.73205080756888;
    -1.73205080756888,1.73205080756888,-1.73205080756888;
    0,1.73205080756888,-1.73205080756888;
    1.73205080756888,1.73205080756888,-1.73205080756888;
    -1.73205080756888,-1.73205080756888,0;
    0,-1.73205080756888,0;
    1.73205080756888,-1.73205080756888,0;
    -1.73205080756888,0,0;
    0,0,0;
    1.73205080756888,0,0;
    -1.73205080756888,1.73205080756888,0;
    0,1.73205080756888,0;
    1.73205080756888,1.73205080756888,0;
    -1.73205080756888,-1.73205080756888,1.73205080756888;
    0,-1.73205080756888,1.73205080756888;
    1.73205080756888,-1.73205080756888,1.73205080756888;
    -1.73205080756888,0,1.73205080756888;
    0,0,1.73205080756888;
    1.73205080756888,0,1.73205080756888;
    -1.73205080756888,1.73205080756888,1.73205080756888;
    0,1.73205080756888,1.73205080756888;
    1.73205080756888,1.73205080756888,1.73205080756888];
L=1; % Length of beam

% n: # of variables
% r: # of points
n = 3;
r = 3^n;

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
% sigmaN_x2 = zeta_x2*mu_x2;
% muN_x2 = mu_x2*(1-log(mu_x2)+lambda_x2)
muN_x2 = log(mu_x2 / sqrt(1 + sigma_x2^2/mu_x2^2));
sigmaN_x2 = log(1 + sigma_x2^2/mu_x2^2);

% I~Normal
muN_x3 = 10;
cov_x3 = 0.05;
sigmaN_x3 = muN_x3*cov_x3;

%% Training points
n_train = 7;
P = linspace( mu_x1 - 3*sigma_x1, mu_x1 + 3*sigma_x1, n_train);
E = linspace( mu_x2 - 3*sigma_x2, mu_x2 + 3*sigma_x2, n_train);
I = linspace( muN_x3 - 3*sigmaN_x3, muN_x3 + 3*sigmaN_x3, n_train);
[PP, EE, II] = ndgrid(P,E,I);
% Input
inp = [PP(:) EE(:) II(:)];
% Output
Y = inp(:,1)./(inp(:,2).*inp(:,3));

% Input in terms of Standard Random Variates
x1 = icdf('Normal',cdf('Extreme Value', inp(:,1), mu_x1, sigma_x1), 0, 1);
x2 = (log(inp(:,2))-lambda_x2)/zeta_x2;
x3 = ((inp(:,3)-muN_x3)/sigmaN_x3);

%% Linear regression analysis to find the trend
X = [ones(size(Y)) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3];
b = regress(Y, X);

%% MCS on the original model
%--------------------------------------------------------------------------
N = 20000;
% Uniform numbers
r = rand(N,3);
x1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1);
x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2);
x3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3);
MCS_disp = x1*L./(x2.*x3);
%-----------------------------------------------------------------------------------------------------
% End of MCS

% Surrogate model - PCE
x1 = icdf('Normal',cdf('Extreme Value', x1, mu_x1, sigma_x1), 0, 1);
x2 = (log(x2)-lambda_x2)/zeta_x2;
x3 = ((x3-muN_x3)/sigmaN_x3);

aa = [ones(N,1) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3];
Y_pce = sum(b'.*[ones(N,1) x1 x2 x3 (x1.^2-1) (x2.^2-1) (x3.^2-1) x1.*x2 x1.*x3 x2.*x3],2);

% # of unknown coeff. for a p order Hermite PCE n random variables
p = 2; n = 3;
Nc = factorial(n+p)/(factorial(n)*factorial(p));





%% Plots %%
%-------------------
% Ploting PDF
%-------------------
figure(1);
% PCE pdf
[fp,xp]=ksdensity(Y_pce);
plot(xp,fp,'-*b')
hold on;
% MCS pdf
[fm,xm] = ksdensity(MCS_disp);
plot(xm,fm,'.-r')
xlabel('Displacement');
ylabel('PDF');
legend('PCE','MCS');
% axis([0 0.3 0 35])
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'PDF', 'pdf') %Save figure


%% Plots %%
%-------------------
% Ploting Histogram
%-------------------
figure(2);
% PCE pdf
subplot(1,2,1);
histfit(Y_pce);
title('Surrogate Model');
xlabel('Displacement');
ylabel('Frequency');
% hold all
% MCS pdf
subplot(1,2,2);
histfit(MCS_disp);
title('Original Model');
xlabel('Displacement');
ylabel('Frequency');
% legend('PCE','MCS');
% axis([0 0.3 0 35])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'Hist','-dpdf'); %Save figure

% %-------------------
% % Ploting CDF
% %-------------------
% figure(3);
% h(1,1) = cdfplot(disp);%PCE
% hold all;
% h(1,2) = cdfplot(MCS_disp);%MCS
% % set( h(1,1), 'LineStyle', '-*', 'Color', 'b');
% % set( h(1,2), 'LineStyle', '.-', 'Color', 'r');
% xlabel('Displacement');
% ylabel('CDF');
% legend('PCE','MCS');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'CDF','-dpdf'); %Save figure
