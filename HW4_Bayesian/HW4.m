% METROPOLIS-HASTINGS BAYESIAN POSTERIOR
rand('seed',12345)
 
% PRIOR OVER SCALE PARAMETERS
B = 1;
 
% DEFINE LIKELIHOOD
likelihood = inline('(B.^A/gamma(A)).*y.^(A-1).*exp(-(B.*y))','y','A','B');
 
% CALCULATE AND VISUALIZE THE LIKELIHOOD SURFACE
yy = linspace(0,10,100);
AA = linspace(0.1,5,100);
likeSurf = zeros(numel(yy),numel(AA));
for iA = 1:numel(AA); likeSurf(:,iA)=likelihood(yy(:),AA(iA),B); end;
 
figure;
surf(likeSurf); ylabel('p(y|A)'); xlabel('A'); colormap hot
 
% DISPLAY CONDITIONAL AT A = 2
hold on; ly = plot3(ones(1,numel(AA))*40,1:100,likeSurf(:,40),'g','linewidth',3)
xlim([0 100]); ylim([0 100]);  axis normal
set(gca,'XTick',[0,100]); set(gca,'XTickLabel',[0 5]);
set(gca,'YTick',[0,100]); set(gca,'YTickLabel',[0 10]);
view(65,25)
legend(ly,'p(y|A=2)','Location','Northeast');
hold off;
title('p(y|A)');
 
% DEFINE PRIOR OVER SHAPE PARAMETERS
prior = inline('sin(pi*A).^2','A');
 
% DEFINE THE POSTERIOR
p = inline('(B.^A/gamma(A)).*y.^(A-1).*exp(-(B.*y)).*sin(pi*A).^2','y','A','B');
 
% CALCULATE AND DISPLAY THE POSTERIOR SURFACE
postSurf = zeros(size(likeSurf));
for iA = 1:numel(AA); postSurf(:,iA)=p(yy(:),AA(iA),B); end;
 
figure
surf(postSurf); ylabel('y'); xlabel('A'); colormap hot
 
% DISPLAY THE PRIOR
hold on; pA = plot3(1:100,ones(1,numel(AA))*100,prior(AA),'b','linewidth',3)
 
% SAMPLE FROM p(A | y = 1.5)
y = 1.5;
target = postSurf(16,:);
 
% DISPLAY POSTERIOR
psA = plot3(1:100, ones(1,numel(AA))*16,postSurf(16,:),'m','linewidth',3)
xlim([0 100]); ylim([0 100]);  axis normal
set(gca,'XTick',[0,100]); set(gca,'XTickLabel',[0 5]);
set(gca,'YTick',[0,100]); set(gca,'YTickLabel',[0 10]);
view(65,25)
legend([pA,psA],{'p(A)','p(A|y = 1.5)'},'Location','Northeast');
hold off
title('p(A|y)');
 
% INITIALIZE THE METROPOLIS-HASTINGS SAMPLER
% DEFINE PROPOSAL DENSITY
q = inline('exppdf(x,mu)','x','mu');
 
% MEAN FOR PROPOSAL DENSITY
mu = 5;
 
% DISPLAY TARGET AND PROPOSAL
figure; hold on;
th = plot(AA,target,'m','Linewidth',2);
qh = plot(AA,q(AA,mu),'k','Linewidth',2)
legend([th,qh],{'Target, p(A)','Proposal, q(A)'});
xlabel('A');
 
% SOME CONSTANTS
nSamples = 5000;
burnIn = 500;
minn = 0.1; maxx = 5;
 
% INTIIALZE SAMPLER
x = zeros(1 ,nSamples);
x(1) = mu;
t = 1;
 
% RUN METROPOLIS-HASTINGS SAMPLER
while t < nSamples
    t = t+1;
 
    % SAMPLE FROM PROPOSAL
    xStar = exprnd(mu);
 
    % CORRECTION FACTOR
    c = q(x(t-1),mu)/q(xStar,mu);
 
    % CALCULATE THE (CORRECTED) ACCEPTANCE RATIO
    alpha = min([1, p(y,xStar,B)/p(y,x(t-1),B)*c]);
 
    % ACCEPT OR REJECT?
    u = rand;
    if u < alpha
        x(t) = xStar;
    else
        x(t) = x(t-1);
    end
end
 
% DISPLAY MARKOV CHAIN
figure;
subplot(211);
stairs(x(1:t),1:t, 'k');
hold on;
hb = plot([0 maxx/2],[burnIn burnIn],'g--','Linewidth',2)
ylabel('t'); xlabel('samples, A');
set(gca , 'YDir', 'reverse');
ylim([0 t])
axis tight;
xlim([0 maxx]);
title('Markov Chain Path');
legend(hb,'Burnin');
 
% DISPLAY SAMPLES
subplot(212);
nBins = 100;
sampleBins = linspace(minn,maxx,nBins);
counts = hist(x(burnIn:end), sampleBins);
bar(sampleBins, counts/sum(counts), 'k');
xlabel('samples, A' ); ylabel( 'p(A | y)' );
title('Samples');
xlim([0 10])
 
% OVERLAY TARGET DISTRIBUTION
hold on;
plot(AA, target/sum(target) , 'm-', 'LineWidth', 2);
legend('Sampled Distribution',sprintf('Target Posterior'))
axis tight