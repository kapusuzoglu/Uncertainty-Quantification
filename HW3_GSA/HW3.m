% Sobol Sensitivity Analysis - GSA
% Author: Berkcan Kapusuzoglu, Vanderbilt University
% 5 Feb 19

% We choose a trial function: Y = x1/(x2*x3) - beam deflection
clear all; close all; clc;
% N=30; % total no. of simulated model output to quantify sensitivity
n_p = 3; % number of model parameters to perform sensitivity analysis on
x=ones(n_p,1);
% n_train = ceil(N/n_p); % base no. of simulations so total no. of model eval ~ N
n_train = 500;
%----------------------------------------
% determine parameter space for computations
% random sampling of parameter space using uniform distribution

% P~Type I
% Type I Extreme Value Distribution: 
mu_x1 = 2000;
cov_x1 = 0.2;
sigma_x1 = mu_x1*cov_x1;

% E~Lognormal
mu_x2 = 30000;
cov_x2 = 0.1;
sigma_x2 = mu_x2*cov_x2;
zeta_E = sqrt(log(1+cov_x2^2));
lam_E = log(mu_x2) - 0.5*zeta_E^2;

% I~Normal
muN_x3 = 10;
cov_x3 = 0.05;
sigmaN_x3 = muN_x3*cov_x3;

dev = 3;
% P = linspace( mu_x1 - dev*sigma_x1, mu_x1 + dev*sigma_x1, n_train);
% E = linspace( mu_x2 - dev*sigma_x2, mu_x2 + dev*sigma_x2, n_train);
% I = linspace( muN_x3 - dev*sigmaN_x3, muN_x3 + dev*sigmaN_x3, n_train);

% PS=[((1 + 1000.*rand(n_train,1))) ((10 + 90.*rand(n_train,1))) ...
%     ((1 + 10.*rand(n_train,1)))];
% 
% % complementary parameter space--need to be defined for MCS
% comp_PS=[((1 + 1000.*rand(n_train,1))) ((10 + 90.*rand(n_train,1))) ...
%     ((1 + 10.*rand(n_train,1)))];

a1 = mu_x1 - dev*sigma_x1;
b1 = mu_x1 + dev*sigma_x1;
a2 = mu_x2 - dev*sigma_x2;
b2 = mu_x2 + dev*sigma_x2;
a3 = muN_x3 - dev*sigmaN_x3;
b3 = muN_x3 + dev*sigmaN_x3;

PS=[(a1 + (b1-a1).*rand(n_train,1)) (a2 + (b2-a2).*rand(n_train,1)) ...
    (a3 + (b3-a3).*rand(n_train,1))];

% complementary parameter space--need to be defined for MCS
comp_PS=[(a1 + (b1-a1).*rand(n_train,1)) (a2 + (b2-a2).*rand(n_train,1)) ...
    (a3 + (b3-a3).*rand(n_train,1))];


output = ones(n_train,1);
c_out_1 = ones(n_train,1);
c_out_t = ones(n_train,1);
for i=1:n_train
    x=PS(i,:);
    output(i,:)=Sobol_obj(x);
    for j=1:n_p
        x = [comp_PS(i,1:j-1),PS(i,j),comp_PS(i,j+1:n_p)];
        c_out_1(i,j) = Sobol_obj(x);
        x = [PS(i,1:j-1),comp_PS(i,j),PS(i,j+1:n_p)];
        c_out_t(i,j) = Sobol_obj(x);
    end
end

% compute variances
n_out = size(output,2); %size of output of objective function
f0 = zeros(1,n_out); % integral of model output
D = zeros(1,n_out); % total variance
% monte carlo integrations to estimate integral functions
for i = 1:n_train
    f0 = f0 + output(i,:)/n_train; % estiamte integral of model output
    D = D + output(i,:).^2/n_train; % start computation of total variance of model output
end
D=D-f0.^2;
Dj=ones(n_p,1)*D; % partial variances associated with parameter j
Dtotj=zeros(n_p,n_out); % total partial variance associated with parameter j


for i = 1:n_train
    for j = 1:n_p
        % start computation of partial variances
        Dj(j,:) = Dj(j,:)-(output(i,:)-c_out_1(i,j)).^2/(2*n_train); 
        % total variance due to pj
        Dtotj(j,:) = Dtotj(j,:)+(output(i,:)-c_out_t(i,j)).^2/(2*n_train);
    end
end
%compute sensitivity indices from variances
Sob_1 = Dj./(ones(n_p,1)*D) %first order
Sob_t = Dtotj./(ones(n_p,1)*D) % total effect



% %sort sensitivity rankings
% %rank_Sob_1j: sorted Sobol first order rankings in terms of ascending magnitude
% %rank_Sob_1jp: parameters ranked in assending order via Sobol first order values
% [rank_Sob_1 rank_Sob_1p] = sort(Sob_1);
% %rank_Sob_tj: sorted Sobol total effect rankings in terms of ascending magnitude
% %rank_Sob_tjp: parameters ranked in assending order via Sobol total effect values
% [rank_Sob_t rank_Sob_tp] = sort(Sob_t);
% %immediately report back
% sens_ind = Sob_1 %report back first order effect sensitivity indices
% % rank_sens = rank_Sob_1 %report back ranked sensitivity indices
% % rank_sens_p = rank_Sob_1p %report back ranked parameters

par_1 = {
'P'            [Sob_1(1)]
'E'            [Sob_1(2)]
'I'            [Sob_1(3)]};

figure()
bar([par_1{:,2}])
set(gca,'XtickLabel',par_1(:,1))
ylabel('First order Effect Indices');
ylim([0 1]);


par_t = {
'P'            [Sob_t(1)]
'E'            [Sob_t(2)]
'I'            [Sob_t(3)]};

figure()
bar([par_t{:,2}])
set(gca,'XtickLabel',par_t(:,1))
ylabel('Total Effect Indices');
ylim([0 1]);


%% Model
function Y = Sobol_obj(x)
% Beam deflection model
x1=x(1); x2=x(2); x3=x(3);
Y =x1/(x2*x3);
% Y= x1+5*x2*x3^2;
end