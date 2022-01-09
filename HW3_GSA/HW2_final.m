% THE CANTILEVER BEAM SUBJECTED TO CONCENTRATED LOAD
%% METHODS: POLYNOMIAL CHAOS, MONTE CARLO SIMULATION & GAUSSIAN PROCESS MODEL
% x1=P, x2=E and x3=I are random fields
% P ~ Type I Extreme Value(m=2000, cov=0.2)
% E ~ Lognormal(m=30000, cov=0.1)
% I ~ Normal(m=10, cov=0.05)
% Create a surrogate model for Y = P/EI using Gaussian Process

clc;
clear all;close all;
tic
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
lambda_x2 = log(mu_x2 / sqrt(1 + sigma_x2^2/mu_x2^2));
zeta_x2 = log(1 + sigma_x2^2/mu_x2^2);
sigmaN_x2 = zeta_E*mu_x2;
muN_x2 = mu_x2*(1-log(mu_x2)+lam_E);

% I~Normal
muN_x3 = 10;
cov_x3 = 0.05;
sigmaN_x3 = muN_x3*cov_x3;

% deviation from the mean
dev1 = 1.15;
dev2 = 1.25;
dev3 = 2.85;

%% Training points
n_train = 7;
P = linspace( mu_x1 - dev1*sigma_x1, mu_x1 + 1.5*sigma_x1, n_train);
E = linspace( mu_x2 - 1.05*sigma_x2, mu_x2 + 1.25*sigma_x2, n_train);
I = linspace( muN_x3 - 1.5*sigmaN_x3, muN_x3 + dev3*sigmaN_x3, n_train);
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
X = [ones(size(Y)) x1 x2 x3];
b = regress(Y, X);
% Residual (error)
residual = Y-X*b;

% x1 = (x1 - min(x1))/(max(x1)-min(x1));
% x2 = (x2 - min(x2))/(max(x2)-min(x2));
% x3 = (x3 - min(x3))/(max(x3)-min(x3));

max(x1)-min(x1)
max(x2)-min(x2)
max(x3)-min(x3)
xi = [x1 x2 x3];

% constant basis (default)
gprMdl = fitrgp(xi,Y,'KernelFunction','ardsquaredexponential');
gprMd2 = fitrgp(xi,Y,'KernelFunction','ardexponential');
% 
% sigma0 = 0.2;
% kparams0 = [3.5, 6.2];
% gprMdl = fitrgp(X,Y,'KernelFunction','squaredexponential',...
%      'KernelParameters',kparams0,'Sigma',sigma0);


%--------------------------------------------------------------------------
%% ----------------------    Gaussian Process     -------------------------

% sigma_e = std(residual);            % std of error
sigma_y = std(Y)/sqrt(2);           % default scalar initial value for the
                                    % noise standard deviation in GP model.
                                    
sigma_x = std(xi)';                 % default length scale vector
% x_initial = [sigma_x; sigma_y];
x_initial = [5 200 .25 0.3];

f = @(x)gpobjective(x,xi,Y);
[xOpt, fval] = fmincon(f,x_initial);

x=xOpt
k_tt = eye(length(Y));
k_tp = zeros(1,length(Y));
for i=1:length(Y)
    for j=1:length(Y)
        k_tt(i,j)= (-0.5*( ((xi(i,1)-xi(j,1))/x(1))^2 + ...
            ((xi(i,2)-xi(j,2))/x(2))^2 + ((xi(i,3)-xi(j,3))/x(3))^2));
    end
end
k_tt = x(4)^2*exp(k_tt);


% MCS on the original model
% -------------------------------------------------------------------------
nsamples = 343;
% Uniform numbers
r = rand(nsamples,3);
p1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1);
e2 = icdf('Lognormal', r(:,2), lambda_x2, zeta_x2);
i3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3);
MCS_disp = p1./(e2.*i3);
% -------------------------------------------------------------------------
% End of MCS

% N_simulations = n_train^3;
y_star = zeros(nsamples,1);
for i=1:nsamples
    pred_pt=randn(1,3);
    for j=1:length(Y)
        k_tp(1,j)= (-0.5*( ((pred_pt(1)-xi(j,1))/x(1))^2 + ...
            ((pred_pt(2)-xi(j,2))/x(2))^2 + ((pred_pt(3)-xi(j,3))/x(3))^2));
    end
    k_tp = x(4)^2*exp(k_tp);
%     k_pp=1;
%     pred_pt=[1 pred_pt];
%     y_star(simulation,1)=k_tp*(k_tt\residual)+pred_pt*b;
    y_star(i,1)=k_tp*(k_tt\Y);
end

% End of GP model
%--------------------------------------------------------------------------


[ypred1, ysd] = resubPredict(gprMdl);
ypred2 = resubPredict(gprMd2);

figure();
histogram(Y)
hold on
histogram(ypred1)
hold on
histogram(y_star)
xlabel('Displacement');
ylabel('Frequency');
legend({'data','MATLAB fitrgp','GPR predictions'},'Location','Best');

figure();
[fy,xy]=ksdensity(Y);
plot(xy,fy,'-r');
hold on
[f1,x1]=ksdensity(ypred1);
plot(x1,f1,'--b');
hold on
[f2,x2]=ksdensity(ypred2);
plot(x2,f2,'-+g');
hold on
[fg,xg]=ksdensity(y_star);
plot(xg,fg,'-*k');
hold on
% MCS pdf
[fm,xm] = ksdensity(MCS_disp);
plot(xm,fm,'.-r')
xlabel('Displacement');
ylabel('PDF');
legend({'Data',['MATLAB fitrgp-' newline 'squared exponential'],...
    ['MATLAB fitrgp-' newline 'exponential'],'GPR predictions','MCS'},'Location','Best');

figure()
plot(Y,'-+r');
hold on
plot(ypred1,'-.b');
hold on
plot(y_star,'-*g');
ylabel('Displacement');
xlabel('# of training points');
legend({'Data',['MATLAB fitrgp-' newline 'squared exponential'],...
    'GPR predictions'},'Location','Best');



toc
