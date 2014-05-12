clc; clear all;

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

disp('LP3 Distribution, Log-Space Moments');
m = mean(log(x));
s = std(log(x));
g = skewness(log(x),0);

% Log-space moment estimators, HH 18.2.27
ahat = 4/g^2;
bhat = 2/(s*g);
zhat = m - 2*s/g;

% Frequency factor Kp, HH 18.2.29
Kp = (2/g)*(1 + g*norminv(0.99)/6 - g^2/36)^3 - 2/g;
q99 = exp(m + s*Kp);

disp(['a = ' num2str(ahat) ', b = ' num2str(bhat) ', z = ' num2str(zhat)]);
disp(['99%: ' num2str(q99) ' m^3/s']);
Kp = (2/g)*(1 + g*norminv((i-3/8)/(n+1/4))/6 - g^2/36).^3 - 2/g;
q = exp(m + s*Kp);
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% 90% KS Bounds (LB Table 7.5)
ca = 0.819/(sqrt(n) - 0.01 + 0.85/sqrt(n));
ub = exp(m + s*((2/g)*(1 + g*norminv((i-1)/n + ca)/6 - g^2/36).^3 - 2/g));
lb = exp(m + s*((2/g)*(1 + g*norminv((i)/n - ca)/6 - g^2/36).^3 - 2/g));
probplot(q,x,lb,ub,'Log-Pearson-3','cms');