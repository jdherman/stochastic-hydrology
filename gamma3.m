clc; clear all;

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

% 3-parameter gamma (LB 7.84)
disp('3-parameter gamma distribution, MM');
that = mean(x) - 2*(sqrt(var(x))/skewness(x,0));
ahat = 4/(skewness(x,0)^2);
bhat = 2/(sqrt(var(x))*skewness(x,0));
disp(['a = ' num2str(ahat) ', b = ' num2str(bhat) ', t = ' num2str(that)]);
disp(['1%: ' num2str(gaminv(.01,ahat,1/bhat)+that) ', 99%: ' num2str(gaminv(.99,ahat,1/bhat)+that) ' m^3/s']);
q = gaminv( (i-3/8)/(n+1/4) , ahat,1/bhat)+that;
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% 90% KS Bounds (LB Table 7.5)
ca = 0.819/(sqrt(n) - 0.01 + 0.85/sqrt(n));
ub = gaminv( (i-1)/n + ca, ahat,1/bhat)+that;
lb = gaminv( (i)/n - ca, ahat,1/bhat)+that;
probplot(q,x,lb,ub,'Gamma-3','cms');