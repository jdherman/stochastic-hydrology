clc; clear all;
addpath mle;

% Gumbel distribution (same as GEV "Type 1" with k=0)
% Note that quantile estimates are linear function of parameters
% so the PPCC are the same regardless of parameter estimates
% Skewness 1.1396

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

% Gumbel MM
disp('Gumbel MM');
ahat = sqrt(var(x))*sqrt(6)/pi;
zhat = mean(x) - 0.5772*ahat;
q1 = zhat - ahat*log(-1*log(0.01));
q99 = zhat - ahat*log(-1*log(0.99));
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat - ahat*log(-1*log((i-3/8)/(n+1/4)));
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% Gumbel LM
disp('Gumbel L-Moments');
xs = sort(x, 'ascend');
suml = 0;
for j=2:n
    suml = suml + (j-1)*xs(j);
end

l2hat = 2/(n*(n-1)) * suml - mean(x);
ahat = l2hat/log(2);
zhat = mean(x) - 0.5772*ahat;
q1 = zhat - ahat*log(-1*log(0.01));
q99 = zhat - ahat*log(-1*log(0.99));
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat - ahat*log(-1*log((i-3/8)/(n+1/4)));
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% Gumbel MLE - need to solve optimization
disp('Gumbel MLE');
options=optimset('Display','off');
[optpt, fval, exitflag] = fsolve(@gumbel_mle, [600,600], options);
ahat = optpt(2);
zhat = optpt(1);
q1 = zhat - ahat*log(-1*log(0.01));
q99 = zhat - ahat*log(-1*log(0.99));
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat - ahat*log(-1*log((i-3/8)/(n+1/4)));
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% Probability plot - note avoid negative logs (complex numbers)
% Normally the inverse functions do this for you but this one is manual
% 90% KS Bounds (LB Table 7.5)
ca = 0.819/(sqrt(n) - 0.01 + 0.85/sqrt(n));
ub = zhat - ahat*log(-1*log((i-1)/n + ca));
lb = zhat - ahat*log(-1*log((i)/n - ca));
lb(imag(lb)~=0) = NaN;
ub(imag(ub)~=0) = NaN;
probplot(q,x,lb,ub,'Gumbel','cms');


