%% DO M-H

% pri_c =@(c,s) 1./(xPrior(:,2)-xPrior(:,1))';
% L =@(c,s) 1/sqrt(2*pi*s)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% posterior dist.

proprnd = @(th) [200+(1000 - 200)*rand, 0.1+(0.5-0.1)*rand];   % proposal random sampler
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);

fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );

f  = @(i,fun) @(th) fun(th, i); %anonymous function whose output is an anonymous function
Likelihood= 1;
for i=1:10
    l1=Likelihood;
    Likelihood =f(i, fun);
    LL = @(th) l1(th).*Likelihood(th);
end

post =@(th) LL(th)*pri_c*pri_s;
% logpost =@(th) log(LL(th))+log(pri_c)+log(pri_s);



tic
trace_mh = mhsample([300 0.05], nSamples,'pdf',post,'proprnd',proprnd, ...
                    'symmetric',true);
toc
%% DO M-H

% pri_c =@(c,s) 1./(xPrior(:,2)-xPrior(:,1))';
% L =@(c,s) 1/sqrt(2*pi*s)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% posterior dist.

proprnd = @(th) [200+(1000 - 200)*rand, 0.1+(0.5-0.1)*rand];   % proposal random sampler
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);

fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );

f  = @(i,fun) @(th) fun(th, i); %anonymous function whose output is an anonymous function
Likelihood= 1;
for i=1:10
    l1=Likelihood;
    Likelihood =f(i, fun);
    LL = @(th) l1(th).*Likelihood(th);
end

post =@(th) LL(th)*pri_c*pri_s;
% logpost =@(th) log(LL(th))+log(pri_c)+log(pri_s);



tic
trace_mh = mhsample([300 0.05], nSamples,'pdf',post,'proprnd',proprnd, ...
                    'symmetric',true);
toc


%% DO M-H

% pri_c =@(c,s) 1./(xPrior(:,2)-xPrior(:,1))';
% L =@(c,s) 1/sqrt(2*pi*s)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% posterior dist.

proprnd = @(th) [200+(1000 - 200)*rand, 0.1+(0.5-0.1)*rand];   % proposal random sampler
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);

fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );

Likelihood=@(th) 1;
for i=1:10
    Likelihood = Likelihood(th).*fun(th,i);
end

post =@(th) Likelihood(th)*pri_c*pri_s;
% logpost =@(th) log(LL(th))+log(pri_c)+log(pri_s);



tic
trace_mh = mhsample([300 0.05], nSamples,'pdf',post,'proprnd',proprnd, ...
                    'symmetric',true);
toc




% pri_c =@(c,s) 1./(xPrior(:,2)-xPrior(:,1))';
% L =@(c,s) 1/sqrt(2*pi*s)^N * exp(-0.5*((yobs-Y(c))*(yobs-Y(c))')/s^2 );
% posterior dist.










proprnd = @(th) [200+(1000 - 200)*rand, 0.1+(0.5-0.1)*rand];   % proposal random sampler
% logpost =@(th) log(L(th))+log(pri_c)+log(pri_s);

fun=@(th,i) 1/sqrt(2*pi*th(2)^2) * exp(-0.5*((yobs(i)-Y(th(1)))^2)/th(2)^2 );

f  = @(i,fun) @(th) fun(th, i); %anonymous function whose output is an anonymous function
Likelihood= 1;
for i=1:10
    l1=Likelihood;
    Likelihood =f(i, fun);
    LL = @(th) l1(th).*Likelihood(th);
end

post =@(th) LL(th)*pri_c*pri_s;
% logpost =@(th) log(LL(th))+log(pri_c)+log(pri_s);



tic
trace_mh = mhsample([300 0.05], nSamples,'pdf',post,'proprnd',proprnd, ...
                    'symmetric',true);
toc
