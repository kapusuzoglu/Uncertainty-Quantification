% THE CANTILEVER BEAM SUBJECTED TO CONCENTRATED LOAD
%% METHODS: POLYNOMIAL CHAOS, MONTE CARLO SIMULATION & GAUSSIAN PROCESS MODEL
% P, E and I are random fields
clc;
clear all;close all;
%---------------------------------------------------------------------------------
%% Polynomial Chaos Expansion Model
%---------------------------------------------------------------------------------
%collocation points
P=[0 0;1.732 0;0 1.732;-1.732 0;0 -1.732;1.732 -1.732;
-1.732 1.732;1.732 1.732;-1.732 -1.732];
L=1; % Length of beam
x=0.25:0.25:1; % Dividing the beam in 4 regions
% --------------------------------------------------------------------------------
%Auto-correlation matrix for EI
for i=1:1:4
    for j=1:1:4
        tau=(j-i)*0.25;
        Cor_ei(i,j)=exp(-abs(tau));
    end
end
Cor_ei = corrcov(Cor_ei);
% KL expansion for EI
[f1,lambda1]=eigs(Cor_ei,2);
for j=1:9
    EI(:,j)=1+0.2*(sqrt(lambda1(1,1)).*f1(:,1)*P(j,1)+sqrt(lambda1(2,2)).*f1(:,2)*P(j,2));
end
%Auto-correlation matrix for w
for i=1:1:4
    for j=1:1:4
        tau=(j-i)*0.25;
        Cor_w(i,j)=exp(-abs(tau));
    end
end
% KL expansion for w
[f2,lambda2]=eigs(Cor_w,2);
for j=1:9
    w(:,j)=0.5+0.1*(sqrt(lambda2(1,1)).*f2(:,1)*P(j,1)+sqrt(lambda2(2,2)).*f2(:,2)*P(j,2));
end
% Analytical expression for tip displacement
for m=1:9
    for n=1:9
        del_tip(m,n)= w(1,n)*L^4/2048*EI(1,m)+...
            (19*w(1,n)/3072+7*w(2,n)/6144)*L^4/EI(2,m)+...
            (61*w(1,n)/3072+31*w(2,n)/3072+11*w(3,n)/6144)*L^4/EI(3,m)+...
            (127*w(1,n)/3072+85*w(2,n)/3072+43*w(3,n)/3072+5*w(4,n)/2048)*L^4/EI(4,m);
    end
end
Y=[del_tip(1,:) del_tip(2,:) del_tip(3,:) del_tip(4,:) del_tip(5,:) del_tip(6,:) del_tip(7,:) del_tip(8,:) del_tip(9,:)]';
% 81 training points
collo1=-1.732;
collo2=0;
collo3=1.732;
xi=PtPoints(collo1,collo2,collo3);
%------------------------------------------------------------------------------
% PCE Surrogate model
for i=1:81
    X(i,:)=[1 xi(i,1) xi(i,2) xi(i,3) xi(i,4) (xi(i,1))^2-1 (xi(i,2))^2-1 (xi(i,3))^2-1 (xi(i,4))^2-1 xi(i,1)*xi(i,2) xi(i,1)*xi(i,3) xi(i,1)*xi(i,4) xi(i,2)*xi(i,3) xi(i,2)*xi(i,4) xi(i,3)*xi(i,4)];
end
b=inv([X'*X])*[X'*Y]; % PCE Surrogate model parameters
% Generate pdf and cdf for PCE model
N_simulations=10000;
for sim=1:N_simulations
    si=randn(1,4);
    tip(:,sim)=b(1)+b(2)*si(1)+b(3)*si(2)+b(4)*si(3)+b(5)*si(4)+...
        b(6)*(1-si(1)^2)+b(7)*(1-si(2)^2)+b(8)*(1-si(3)^2)+b(9)*(1-si(4)^2)+...
        b(10)*si(1)*si(2)+b(11)*si(1)*si(3)+b(12)*si(1)*si(4)+b(13)*si(2)*si(3)+b(14)*si(2)*si(4)+b(15)*si(3)*si(4);
end
%
%--------------------------------------------------------------------------------------------------------------
% End of PCE


%% MCS on the original model
%--------------------------------------------------------------------------------------------------------------
for m=1:200
    for n=1:200
        w=0.5+0.1*(sqrt(lambda2(1,1))*f2(:,1)*randn+sqrt(lambda2(2,2))*f2(:,2)*randn);
        EI=1+0.2*(sqrt(lambda1(1,1))*f1(:,1)*randn+sqrt(lambda1(2,2))*f1(:,2)*randn);
        mcs_tip(m,n)= w(1)*L^4/2048*EI(1)+...
            (19*w(1)/3072+7*w(2)/6144)*L^4/EI(2)+...
            (61*w(1)/3072+31*w(2)/3072+11*w(3)/6144)*L^4/EI(3)+...
            (127*w(1)/3072+85*w(2)/3072+43*w(3)/3072+5*w(4)/2048)*L^4/EI(4);
    end
end
Z=[mcs_tip(:)]';
%-----------------------------------------------------------------------------------------------------
% End of MCS


%% Gaussian Process
%-----------------------------------------------------------------------------------------------------
x_initial=[1.2 1.4 1.1 1.5 0.3];
f = @(x)gpobjective(x,Y);
[xOpt, fval] = fminunc(f,x_initial);
%Regression Analysis to find TREND
X_trend=[ones(81,1) xi];
b_trend=(X_trend'*X_trend)\(X_trend'*Y);
%Residaul
residual=Y-X_trend*b_trend;
x=xOpt;
for i=1:81
    for j=1:81
        k_tt(i,j)= x(5)^2*exp(-0.5*((xi(i,1)-xi(j,1))^2/x(1)^2+(xi(i,2)-xi(j,2))^2/x(2)^2+...
            (xi(i,3)-xi(j,3))^2/x(3)^2+(xi(i,4)-xi(j,4))^2/x(4)^2));
    end
end
for simulation=1:N_simulations
    pred_pt=randn(1,4);
    for j=1:81
        k_tp(1,j)= x(5)^2*exp(-0.5*((pred_pt(1)-xi(j,1))^2/x(1)^2+(pred_pt(2)-xi(j,2))^2/x(2)^2+...
            (pred_pt(3)-xi(j,3))^2/x(3)^2+(pred_pt(4)-xi(j,4))^2/x(4)^2));
    end
    k_pp=1;
    pred_pt=[1 pred_pt];
    y_star(simulation,1)=k_tp*(inv(k_tt)*residual)+pred_pt*b_trend;
end
%-----------------------------------
% End of GP model

%% Plots %%
%-------------------
% Ploting PDF
%-------------------
figure(1);
% PCE pdf
%hist(tip)
[fp,xp]=ksdensity(tip);
plot(xp,fp,'--r')
hold on;
% MCS pdf
%hist(MCS_tip)
[fm,xm] = ksdensity(Z);
plot(xm,fm,'--b')
hold on;
% GP pdf
[fg,xg]=ksdensity(y_star);
plot(xg,fg,'-*k');
hold on;
xlabel('Tip displacement');
ylabel('PDF');
legend('PCE','MCS');
axis([0 0.3 0 35])

% %-------------------
% % Ploting CDF
% %-------------------
% figure(2);
% cdfplot(tip);%PCE
% hold on;
% cdfplot(Z);%MCS
% hold on;
% cdfplot(y_star);%GP
% hold on;
% xlabel('Tip displacement');
% ylabel('CDF');
% legend('PCE','MCS','GP');
