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
zeta_E = sqrt(log(1+cov_x2^2));
lam_E = log(mu_x2) - 0.5*zeta_E^2;
sigmaN_x2 = zeta_E*mu_x2;
muN_x2 = mu_x2*(1-log(mu_x2)+lam_E);

% I~Normal
muN_x3 = 10;
cov_x3 = 0.05;
sigmaN_x3 = muN_x3*cov_x3;


% -------------------------------------------------------------------------
% Auto-correlation matrix for x1, x2, x3
N=3;
tau = linspace(0,1,N);
rxx =  exp(-abs(-tau));
Cor_x1 = toeplitz(rxx);
% KL expansion for x1, x2, x3
[f1,lambda1]=eigs(Cor_x1);
f2 = f1; f3 = f1; lambda2 = lambda1; lambda3 = lambda1;
for j=1:r
    x1(:,j) = muN_x1 + (sqrt(lambda1(1,1)).*f1(:,1)*P(j,1)+...
        sqrt(lambda1(2,2)).*f1(:,2)*P(j,2)+sqrt(lambda1(3,3)).*f1(:,3)*P(j,3));
    x2(:,j) = muN_x2 + (sqrt(lambda1(1,1)).*f2(:,1)*P(j,1)+...
        sqrt(lambda2(2,2)).*f2(:,2)*P(j,2)+sqrt(lambda2(3,3)).*f2(:,3)*P(j,3));
    x3(:,j) = muN_x3 + (sqrt(lambda3(1,1)).*f3(:,1)*P(j,1)+...
        sqrt(lambda3(2,2)).*f3(:,2)*P(j,2)+sqrt(lambda3(3,3)).*f3(:,3)*P(j,3));
end

% Analytical expression for tip displacement
% for i=1:n
%     for j=1:n
%         for k=1:n
%             del_tip(i,j,k) = x1(1,i)*L/(x2(1,j)*x3(1,k));
%         end
%     end
% end
% Y=[del_tip(1,:) del_tip(2,:) del_tip(3,:)]';

del_tip = x1*L./(x2.*x3);
Y = del_tip';

% # of unknown coeff. for a p order Hermite PCE n random variables
p = 2; n = 3;
Nc = factorial(n+p)/(factorial(n)*factorial(p));

% r = 27 training points
collo1=-sqrt(3);
collo2=0;
collo3=sqrt(3);
xi=TrPoints(collo1,collo2,collo3);

    
%------------------------------------------------------------------------
% PCE Surrogate model
X = zeros(r,Nc);
X(:,:) = [ones(r,1) xi(1:r,1) xi(1:r,2) xi(1:r,3) (xi(1:r,1)).^2-1 (xi(1:r,2)).^2-1 ...
    (xi(1:r,3)).^2-1 xi(1:r,1).*xi(1:r,2) xi(1:r,1).*xi(1:r,3) xi(1:r,2).*xi(1:r,3)];
% PCE Surrogate model parameters
b = (X'*X)\(X'*Y);


% Generate pdf and cdf for PCE model
N_simulations = 5000;
psi = randn(N_simulations,n);
psi = [ones(N_simulations,1) psi(:,1) psi(:,2) psi(:,3) (1-psi(:,1).^2) (1-psi(:,2).^2) ...
           (1-psi(:,3).^2) psi(:,1).*psi(:,2) psi(:,1).*psi(:,3) psi(:,2).*psi(:,3)];
    % disp = zeros(N_simulations,1);
disp = b'*psi';
% disp_tot = sum(disp);

% % b=mean(b')';
% for i=1:N_simulations
%     si = randn(1,n);
%     tip(:,i) = b(1)+b(2)*si(1)+b(3)*si(2)+b(4)*si(3)+ b(5)*(1-si(1)^2)+...
%         b(6)*(1-si(2)^2)+b(7)*(1-si(3)^2)+b(8)*si(1)*si(2)+b(9)*si(1)...
%         *si(3)+b(10)*si(2)*si(3);
% end
%
%--------------------------------------------------------------------------
% End of PCE
% 
% 
% %% MCS on the original model
% %--------------------------------------------------------------------------
% nsamples = 5000;
% x1 = -log(log(1./rand(nsamples,1)))./alphan + u;
% x2 = exp(lam_E + zeta_E*norminv(rand(nsamples,1)));
% x3 = muN_x3 + sigmaN_x3*norminv(rand(nsamples,1));
% MCS_disp = x1.*L./(x2.*x3);
% Z = MCS_disp(:)';
% %-----------------------------------------------------------------------------------------------------

%% MCS on the original model
%--------------------------------------------------------------------------
nsamples = 10000;
% Uniform numbers
muN_x2 = log(mu_x2 / sqrt(1 + sigma_x2^2/mu_x2^2));
sigmaN_x2 = log(1 + sigma_x2^2/mu_x2^2);
r = rand(nsamples,3);
x1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1);
x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2);
x3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3);
MCS_disp = x1*L./(x2.*x3);
%-----------------------------------------------------------------------------------------------------
% End of MCS


%% Plots %%
%-------------------
% Ploting PDF
%-------------------
figure(1);
% PCE pdf
disp=disp(:);
[fp,xp]=ksdensity(disp);
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
histfit(disp);
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

%-------------------
% Ploting CDF
%-------------------
figure(3);
h(1,1) = cdfplot(disp);%PCE
hold all;
h(1,2) = cdfplot(MCS_disp);%MCS
% set( h(1,1), 'LineStyle', '-*', 'Color', 'b');
% set( h(1,2), 'LineStyle', '.-', 'Color', 'r');
xlabel('Displacement');
ylabel('CDF');
legend('PCE','MCS');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'CDF','-dpdf'); %Save figure
