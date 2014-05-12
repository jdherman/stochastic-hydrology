clc; clear all;

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

% 2-parameter gamma (LB 7.85)
disp('2-parameter gamma distribution, MM');
ahat = mean(x)^2/var(x);
bhat = mean(x)/var(x);
disp(['a = ' num2str(ahat) ', b = ' num2str(bhat)]);
disp(['1%: ' num2str(gaminv(.01,ahat,1/bhat)) ', 99%: ' num2str(gaminv(.99,ahat,1/bhat)) ' m^3/s']);
q = gaminv( (i-3/8)/(n+1/4) , ahat,1/bhat);
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% 90% KS Bounds (LB Table 7.5)
ca = 0.819/(sqrt(n) - 0.01 + 0.85/sqrt(n));
ub = gaminv( (i-1)/n + ca, ahat,1/bhat);
lb = gaminv( (i)/n - ca, ahat,1/bhat);
probplot(q,x,lb,ub,'Gamma-2','cms');