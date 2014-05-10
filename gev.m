clc; clear all;
addpath mle;

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

% GEV LM
disp('GEV L-Moments');
xs = sort(x, 'ascend');

% Calculate b values
b0 = mean(x);
b1 = 0;
for j=2:n
    b1 = b1 + (j-1)*xs(j);
end
b1 = b1/(n*(n-1));
b2 = 0;
for j=3:n
    b2 = b2 + (j-1)*(j-2)*xs(j);
end
b2 = b2/(n*(n-1)*(n-2));

% L-moments and estimators (LB 7.94)
l1 = b0;
l2 = 2*b1 - b0;
l3 = 6*b2 - 6*b1 + b0;
t3 = l3/l2;

c = 2/(3+t3) - log10(2)/log10(3);
khat = 7.8590*c + 2.9554*c^2;
ahat = l2*khat/((1-2^(-khat))*gamma(1+khat));
zhat = l1 - (ahat/khat)*(1 - gamma(1+khat));

q1 = zhat + (ahat/khat)*(1 - (-log(0.01))^khat);
q99 = zhat + (ahat/khat)*(1 - (-log(0.99))^khat);
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat) ', kappa = ' num2str(khat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat + (ahat/khat)*(1 - (-log((i-3/8)/(n+1/4))).^khat);
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% GEV MLE
disp('GEV MLE');
optpt = fminsearch(@gev_mle, [ahat,khat,zhat]);
ahat = optpt(1);
khat = optpt(2);
qhat = optpt(3);

% transform back
if(khat < 0)
    zhat = min(x) - ahat/khat - exp(qhat);
else
    zhat = max(x) - ahat/khat + exp(qhat);
end
    
q1 = zhat + (ahat/khat)*(1 - (-log(0.01))^khat);
q99 = zhat + (ahat/khat)*(1 - (-log(0.99))^khat);
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat) ', kappa = ' num2str(khat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat + (ahat/khat)*(1 - (-log((i-3/8)/(n+1/4))).^khat);
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% GEV GMLE
disp('GEV GMLE');
optpt = fminsearch(@gev_gmle, [ahat,-0.1,zhat]);
ahat = optpt(1);
khat = optpt(2);
qhat = optpt(3);

% transform back
if(khat < 0)
    zhat = min(x) - ahat/khat - exp(qhat);
else
    zhat = max(x) - ahat/khat + exp(qhat);
end

q1 = zhat + (ahat/khat)*(1 - (-log(0.01))^khat);
q99 = zhat + (ahat/khat)*(1 - (-log(0.99))^khat);
disp(['alpha = ' num2str(ahat) ', xi = ' num2str(zhat) ', kappa = ' num2str(khat)]);
disp(['1%: ' num2str(q1) ', 99%: ' num2str(q99) ' m^3/s']);
q = zhat + (ahat/khat)*(1 - (-log((i-3/8)/(n+1/4))).^khat);
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% Probability plot - note avoid negative logs (complex numbers)
% Normally the inverse functions do this for you but this one is manual
ub = zhat + (ahat/khat)*(1 - (-log((i-1)/n + 0.127)).^khat);
lb = zhat + (ahat/khat)*(1 - (-log((i)/n - 0.127)).^khat);
lb(imag(lb)~=0) = NaN;
ub(imag(ub)~=0) = NaN;
probplot(q,x,lb,ub,'GEV','cms');
