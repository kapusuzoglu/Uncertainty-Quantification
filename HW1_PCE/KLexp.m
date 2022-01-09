% HW 1
clc; clear;

t = linspace(0,1,21);
R = exp(-abs(t)); % given autocorrelation function

% autocorrelation function of the random process
Rtt = toeplitz(R);

[V,D] = eig(Rtt);

[d,ind] = sort(diag(D),'descend');
Ds = D(ind,ind); Vs = (V(:,ind)); %sorted

N = 10000;
% SimulatedSignal = zeros(N,1);
for i=1:N
    
   SRVs = normrnd(0,1,[1,11]);
   SimulatedSignal(i,:)=sum(sqrt(transpose(d(1:11))) .* Vs(:,1:11) ...
       .* SRVs, 2);
end

% the 1st row is what we want to compare with
R_sim = corrcoef(SimulatedSignal);

figure()
plot(t,R,'k',t,R_sim(1,:),'k--');
title('KL expansion','interpreter','latex','fontsize',16);
xlabel('y','interpreter','latex','fontsize',16);
ylabel('R(t)','interpreter','latex','fontsize',16);
l1 = legend('Target','Simulated');
set(l1,'interpreter','latex','fontsize',14);

