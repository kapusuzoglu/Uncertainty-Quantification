% clc;
% clear;
close all
% 
% %% Assume uncorrelated inputs
% % Grid Samples
% % P~Type I
% mu_x1 = 2000;
% cov_x1 = 0.2;
% sigma_x1 = mu_x1*cov_x1;
% 
% % E~Lognormal
% mu_x2 = 30000;
% cov_x2 = 0.1;
% sigma_x2 = mu_x2*cov_x2;
% zeta_x2 = sqrt(log(1+cov_x2^2));
% lambda_x2 = log(mu_x2) - 0.5*zeta_x2^2;
% 
% % I~Normal
% muN_x3 = 10;
% cov_x3 = 0.05;
% sigmaN_x3 = muN_x3*cov_x3;
% 
% dev = 3;
% N_train=5;
% P = linspace( mu_x1 - dev*sigma_x1, mu_x1 + dev*sigma_x1, N_train);
% E = linspace( mu_x2 - dev*sigma_x2, mu_x2 + dev*sigma_x2, N_train);
% I = linspace( muN_x3 - dev*sigmaN_x3, muN_x3 + dev*sigmaN_x3, N_train);
% [PP, EE, II] = ndgrid(P,E,I);
% 
% % Input
% inp = [PP(:) EE(:) II(:)];
% 
% % Output
% % Y = inp(:,1)./(inp(:,2).*inp(:,3));
% % 
% %% MCS
% %--------------------------------------------------------------------------
% nsamples = 10;
% % Uniform numbers
% muN_x2 = log(mu_x2 / sqrt(1 + sigma_x2^2/mu_x2^2));
% sigmaN_x2 = log(1 + sigma_x2^2/mu_x2^2);
% r = rand(nsamples,3);
% x1 = icdf('Extreme Value', r(:,1), mu_x1, sigma_x1)';
% x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2)';
% x3 = icdf('Normal', r(:,3), muN_x3, sigmaN_x3)';
% Y =@(c) c.*mu_x1./(mu_x2.*muN_x3);
% %--------------------------------------------------------------------------
% 
% P = mu_x1; E = mu_x2; I = muN_x3;
% L = 10;
% % beam deflection Euler-Bernoulli
% mu_y = P*L^3/(3*E*I);
% sigma_y = 0.33;
% % Generate 10 samples of Y~N(mu_y,sigma_y)
% N = 10;
% yobs = normrnd(mu_y,sigma_y,1,N);
% 
% %% Case I
% % x1, x2, x3 are deterministic at their mean values
% 
% %% Initialize the Metropolis sampler
% T = 10000; % Set the maximum number of iterations
% sigma_c = 35; % Set standard deviation of normal proposal density
% sigma_eps = 0.3; % Set standard deviation of normal proposal density
% cmin = 0; cmax = L^3; % define a range for starting values
% c = zeros( 1 , T ); % Init storage space for our samples
% sigma_e = zeros( 1 , T ); % Init storage space for our samples
% seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
% 
% % Prior distribution for c = L^3/3
% % Assume uniform (non-informative prior)
% c(1) = unifrnd( cmin , cmax ); % Generate start value
% pri_c = 1/(cmax-cmin);
% 
% smin = 0; smax = 1; % define a range for starting values
% sigma_e(1) = unifrnd( smin , smax ); % Generate start value
% % err =@(s) normrnd(0,s,1,N);
% pri_s = 1/(smax-smin);
% 
% % likelihood dist.
% L =@(c,s) 1/sqrt(2*pi*s^2)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% 
% % posterior dist.
% post =@(c,s) L(c,s)*pri_c*pri_s;
% 
% 
% %% Start sampling Metropolis
% t = 1;
% while t < T % Iterate until we have T samples
%     t = t + 1;
%     % Propose a new value for theta using a normal proposal density
%     c_star = normrnd( c(t-1) , sigma_c );
%     sigma_star = normrnd( sigma_e(t-1) , sigma_eps );
%     % Calculate the acceptance ratio
%     alpha = min( [ 1 post( c_star, sigma_star ) / post( c(t-1),sigma_e(t-1) ) ] );
%     % Draw a uniform deviate from [ 0 1 ]
%     u = rand;
%     % Do we accept this proposal?
%     if u < alpha
%         c(t) = c_star; % If so, proposal becomes new state
%         sigma_e(t) = sigma_star; % If so, proposal becomes new state
%     else
%         c(t) = c(t-1); % If not, copy old state
%         sigma_e(t) = sigma_e(t-1); % If so, proposal becomes new state
%     end
% end
% 
% 
% %% Start sampling M-H
% init=[0 0.2];
% 
% % posterior dist.
% logpost =@(c,s) log(L(c,s))+log(pri_c)+log(pri_s);
% 
% D = numel(init);
% samples = zeros(D, T);
% state = init;
% Lp_state = logpost(state(1),state(2));
% for ss = 1:T
%     % Propose
%     prop = state + [sigma_c sigma_eps].*randn(size(state));
%     Lp_prop = logpost(prop(1),prop(2));
%     if log(rand) < (Lp_prop - Lp_state)
%         % Accept
%         state = prop;
%         Lp_state = Lp_prop;
%     end
%     samples(:, ss) = state(:);
% end
% 
% %% uniform prior for c
% X_c = linspace( cmin , cmax , 100 );
% y_c = unifpdf(X_c,cmin,cmax);
% % uniform prior for sigma
% X_s = linspace( smin , smax , 100 );
% y_s = unifpdf(X_s,smin,smax);
% 
% %% DO SLICE SAMPLING
% % likelihood dist.
% L =@(th) 1/sqrt(2*pi*th(2)^2)^N * exp(-0.5*((yobs-Y(th(1)))*(yobs-Y(th(1)))')/th(2)^2 );
% fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );
% 
% f  = @(i,fun) @(th) fun(th, i); %anonymous function whose output is an anonymous function
% Likelihood= 1;
% for i=1:10
%     l1=Likelihood;
%     Likelihood =f(i, fun);
%     LL = @(th) l1(th).*Likelihood(th);
% end
% 
% post =@(th) LL(th)*pri_c*pri_s;
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);
% initial_samples = [300,0.5];
% nSamples = 1e4;
% 
% trace = slicesample(initial_samples, 8000,'logpdf',logpost,'thin',5,...
% 	'burnin', 2000);
% mean(trace(:,2))
% 
% 
% %% DO M-H
% 
% % pri_c =@(c,s) 1./(xPrior(:,2)-xPrior(:,1))';
% % L =@(c,s) 1/sqrt(2*pi*s)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% % posterior dist.
% 
% proprnd = @(th) [200+(1000 - 200)*rand, 0+(0.6-0.1)*rand];   % proposal random sampler
% % logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);
% 
% fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );
% 
% f  = @(i,fun) @(th) fun(th, i); %anonymous function whose output is an anonymous function
% Likelihood= 1;
% for i=1:10
%     l1=Likelihood;
%     Likelihood =f(i, fun);
%     LL = @(th) l1(th).*Likelihood(th);
% end
% 
% post =@(th) L(th)*pri_c*pri_s;
% % logpost =@(th) log(LL(th))+log(pri_c)+log(pri_s);
% 
% % posterior dist.
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);
% 
% unsymPropPdf=@(x,y) mvnpdf([x(1), x(2)],[y(1),y(2)],[0.01*y(1)^2 0;
%                                                          0 0.16*y(2)^2]);
% tic
% trace_mh = mhsample([300 0.1], 1e4,'logpdf',logpost,'proprnd',proprnd, ...
%                     'proppdf',unsymPropPdf);
% toc
% 
% 
% 
% post =@(th) LL(th)*pri_c*pri_s;
% 
% 
% %% Gibbs Sampling
% scale = [10 .01];
% v_th=Gibbs(post,[30 0.4],T,scale);
% 
% 
% % 
% %% Display ksdensity of c
% figure( 1 ); clf;
% subplot( 3,1,1 );
% [f1,c1]=ksdensity(c);
% [f2,c2]=ksdensity(samples(1,:));
% [f3,c3]=ksdensity(trace_mh(:,1));
% [f4,c4]=ksdensity(trace(:,1));
% [f5,c5]=ksdensity(v_th(:,1));
% 
% % set(groot,'DefaultAxesColorOrder',[0 0 0],...
% %       'DefaultAxesLineStyleOrder','-|--|:|-.|-*|o')
% set(groot,'DefaultAxesColorOrder',[0 0 0],...
%       'DefaultAxesLineStyleOrder','-|--|:|-.|-*|-o|-.*')
% plot(X_c,y_c);hold on;
% plot(c1,f1);hold on;
% plot(c2,f2);hold on;
% plot(c3,f3);hold on;
% plot(c4,f4);hold on;
% plot(c5,f5);
% 
% xlim( [ cmin cmax ] );
% xlabel( 'c' ); ylabel( 'PDF' );
% legend('Prior of c','Metropolis','M-H', 'Matlab built-in M-H',...
%     'SliceSampling','Gibbs Sampling')
% 
% 
% %% Display history of our samples
% subplot( 3,1,2:3 );
% stairs( c , 1:T , 'k-' );
% ylabel( '# of iterations' ); xlabel( 'c' );
% set( gca , 'YDir' , 'reverse' );
% xlim( [ cmin cmax ] );
% 
% 
% %% Display ksdensity of sigma
% figure( 2 ); clf;
% subplot( 3,1,1 );
% plot(X_s,y_s);hold on;
% ksdensity(sigma_e);hold on;
% ksdensity(samples(2,:));hold on;
% ksdensity(trace_mh(:,2));hold on;
% ksdensity(trace(:,2));hold on;
% ksdensity(v_th(:,2));
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'PDF' );
% legend('Prior of \sigma_{obs}','Metropolis','M-H','Matlab built-in M-H',...
%    'SliceSampling','Gibbs Sampling')
% 
% 
% subplot( 3,1,2 );
% mu = 0;
% sigma = mean(trace_mh(:,2));
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% x =-2:0.1:2;
% y = pdf(pd,x);
% plot(x,y,'LineWidth',2);
% xlabel( '\epsilon' ); ylabel( 'PDF' );
% 
% 
% subplot( 3,1,3 );
% mu = mean(trace_mh(:,1))*P/E/I;
% sigma = std(trace_mh(:,2));
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% x =1.5:0.01:3;
% y = pdf(pd,x);
% plot(x,y,'LineWidth',2);
% xlabel( 'Y' ); ylabel( 'PDF' );
% 
% figure()
% subplot(3,2,3)
% histogram(sigma_e)
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('Metropolis')
% subplot(3,2,1)
% histogram(samples(2,:));
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('M-H')
% subplot(3,2,2)
% histogram(trace_mh(:,2));
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('Matlab built-in M-H')
% subplot(3,2,4)
% histogram(trace(:,2));
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('SliceSampling')
% subplot(3,2,5)
% histogram(v_th(:,2));
% xlim( [ smin 1 ] );
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('Gibbs Sampling')
% 
% 
% % 
% % %% Display history of our samples
% % subplot( 3,1,2:3 );
% % stairs( sigma_e , 1:T , 'k-' );
% % ylabel( '# of iterations' ); xlabel( '\sigma_{obs}' );
% % set( gca , 'YDir' , 'reverse' );
% % % xlim( [ smin 1 ] );
% 








%% Case II
% x1, x2, x3 have given distributions
% To calculate likelihood assume error_obs~N(0,sigma_obs)
% We have 2 calibration quantity: c and error_obs


% Generate 10 samples of Y~N(mu_y,sigma_y)
N = 10;
yobs = normrnd(2.2,0.33,1,N);

% log posterior
logpost =@(X) ProposedPosterior(X,yobs,0,1000,0,1);
% proposed random number generator
proprnd = @(X) mvnrnd([500,0.4],[150^2 0;0 0.12^2],1);
% proposed unsymmetric pdf
unsymmPropPdf=@(X,Y) mvnpdf([X(1), X(2)],[Y(1),Y(2)],[0.09*Y(1)^2 0;
                                                         0 0.09*Y(2)^2]);

%                                                      
% tic
% trace_mm = mhsample([500 0.4], 1000,'logpdf',logpost,'proprnd',proprnd, ...
%                     'symmetric',true,'thin',5,'burnin', 1000);
% toc
% 
% tic
% trace_mh = mhsample([500 0.4], 1e3,'logpdf',logpost,'proprnd',proprnd, ...
%                     'proppdf',unsymmPropPdf,'thin',5,'burnin', 1000);
% toc
% 
% figure()
% subplot(1,2,1)
% histogram(trace_mm(:,1))
% xlabel( 'c' ); ylabel( 'Frequency' );
% legend('Metropolis')
% subplot(1,2,2)
% histogram(trace_mm(:,2))
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('Metropolis')
% 
% figure()
% subplot(1,2,1)
% histogram(trace_mh(:,1))
% xlabel( 'c' ); ylabel( 'Frequency' );
% legend('M-H')
% subplot(1,2,2)
% histogram(trace_mh(:,2))
% xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
% legend('M-H')


tic
trace = slicesample([500 0.4], 100,'logpdf',logpost,'thin',5,...
	'burnin', 2000,'width',[100 0.4]);
mean(trace(:,2))
toc

figure()
subplot(1,2,1)
histogram(trace(:,1))
xlabel( 'c' ); ylabel( 'Frequency' );
legend('SliceSampling')
subplot(1,2,2)
histogram(trace(:,2))
xlabel( '\sigma_{obs}' ); ylabel( 'Frequency' );
legend('SliceSampling')





function Y = ProposedPosterior(X, y_obs, cmin, cmax, sigma_min, sigma_max)
% c = X(1), sigma_obs = X(2)
N = 500; % # of MCS for X
% Uniform numbers (LHS)
r = [lhsdesign(N,3) rand(N,1)];
u = 0.5772;
sigmaBar = sqrt(6)*(2000*0.2)/pi;
muBar=2000-u*sigmaBar;
muN_x2 = log(30000 / sqrt(1 + 3000^2/30000^2));
sigmaN_x2 = log(1 + 3000^2/30000^2);
%% MCS
%--------------------------------------------------------------------------
x1 = icdf('GeneralizedExtreme Value', r(:,1), sigmaBar, muBar);
x2 = icdf('Lognormal', r(:,2), muN_x2, sigmaN_x2);
x3 = icdf('Normal', r(:,3), 10, 0.5);

% Y_obs = Y_m - eps_obs
eps_obs = icdf('Normal', r(:,4), 0, abs( X(2) ));
Y_m = X(1)*x1./(x2.*x3);
Y_obs = Y_m - eps_obs;
% find pdf values of corresponding y_obs
[fobs, yobs] = ksdensity(Y_obs);
len = length(y_obs);
L=0;
for i=1:len
   [~,idx] = min(abs(yobs-y_obs(i))); 
   pdf_value = fobs(idx);
   L = L + log(pdf_value);
end


Y = L + log(pdf('Uniform',X(1),cmin,cmax))+log(pdf('Uniform',X(2),sigma_min,sigma_max));
end

