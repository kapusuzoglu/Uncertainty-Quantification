%% clear memory, screen, and close all figures
clear, clc, close all;

%% Inputs
% P~Type I
mu_x1 = 2000;
% E~Lognormal
mu_x2 = 30000;
% I~Normal
mu_x3 = 10;

%% Inputs for the Observation/Measurement/Process Model
L = 10; % length of the beam
% mean and std of beam deflection Euler-Bernoulli
mu_y = mu_x1*L^3/(3*mu_x2*mu_x3); sigma_y = 0.3;

% Y_obs = Y_m - eps_obs
y_obs = normrnd(mu_y,sigma_y,1,10);
% # of obs.
obs = length(y_obs);

% # of particles
N = 3000;

% Uniform numbers
u = rand(N,2);

% Initial, prior particles
initParticles = zeros(N,2);
initParticles(:,1) = 1000 * u(:,1);
initParticles(:,2) = u(:,2);

% Initial weight for each particle
w = 1/N * ones(N,1);

Particles = zeros(N,2,obs+1);
Particles(:,:,1) = initParticles;

for i=1:obs
    tic
   CurrState = Particles(:,:,i);
   w = myWeight(y_obs(i), CurrState(:,1), CurrState(:,2));
   Particles(:,:,i+1) = Resampling(CurrState, w);
    toc
end


%% Figures
figure
[p,f] = ksdensity(Particles(:,1,1));
[p1,f1] = ksdensity(Particles(:,1,2));
[p5,f5] = ksdensity(Particles(:,1,5));
[p10,f10] = ksdensity(Particles(:,1,11));
plot(f,p,'-r',f1,p1,'--g',f5,p5,'-.b',f10,p10,'-*k');
title('Evolution of C - Deterministic input','FontSize',14)
xlabel('C');
ylabel('PDF');
legend('Prior','Iteration 1','Iteration 5','Iteration 10')


%% Figures
figure
[p,f] = ksdensity(Particles(:,2,1));
[p1,f1] = ksdensity(Particles(:,2,2));
[p5,f5] = ksdensity(Particles(:,2,5));
[p10,f10] = ksdensity(Particles(:,2,11));
plot(f,p,'-r',f1,p1,'--g',f5,p5,'-.b',f10,p10,'-*k');
title('Evolution of \sigma_{obs} - Deterministic input','FontSize',14)
xlabel('\sigma_{obs}');
ylabel('PDF');
legend('Prior','Iteration 1','Iteration 5','Iteration 10')
xlim([0 1])

function wk = myWeight(y_obs, C, sigma_obs)

Nn = length(C);

w = zeros(Nn,1);

for i=1:Nn
    C_i = C(i,1);
    sigma_obs_i = sigma_obs(i,1);
    
    eps_obs = normrnd(0,sigma_obs_i);
    Y_m = C_i*2000/(30000*10);
    Y_obs = Y_m - eps_obs
    % find pdf values of corresponding y_obs
    [fobs, yobs] = ksdensity(Y_obs);

    [~,ind] = min(abs(yobs-y_obs)); 
    w(i,1) = fobs(ind);

end
wk = w/sum(w);
end


function NextParticles = Resampling(CurrState, w)

% n = length(CurrState);
% edges = min([0 cumsum(w)'],1); % protect against accumulated round-off
% edges(end) = 1;                 % get the upper edge exact
% u1 = rand/n;
% 
% % this works like the inverse of the empirical distribution and returns
% % the interval where the sample is to be found
% [~, idx] = histc(u1:1/n:1, edges);
% Particles = CurrState(idx,:);                   % extract new particles

NextParticles = zeros(size(CurrState));

len = length(w);

u = rand(len,1);
myCumsum = cumsum(w);
myCumsum(end) = 1;

for i=1:len
    [~,ind] = min(abs(myCumsum - u(i)));
    NextParticles(i,:) = CurrState(ind,:);
end

end

