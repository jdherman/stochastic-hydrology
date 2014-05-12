clc; clear all;

x = load('input_data.txt');

% Normal distribution
m = mean(x); v = var(x);
disp('Normal distribution (MM/MLE)');
disp(['m = ' num2str(m) ', s2 = ' num2str(v) ' (s = ' num2str(sqrt(v)) ')']);
disp(['1%: ' num2str(norminv(.01,m,sqrt(v))) ', 99%: ' num2str(norminv(.99,m,sqrt(v))) ' m^3/s']);

% Probability Plot w/90% KS Bounds
i = 1:length(x);
n = length(x);
q = norminv( (i-3/8)/(n+1/4) , m, sqrt(v));
% 90% KS Bounds (LB Table 7.5)
ca = 0.819/(sqrt(n) - 0.01 + 0.85/sqrt(n));
ub = norminv( (i-1)/n + ca, m, sqrt(v));
lb = norminv( (i)/n - ca, m, sqrt(v));
probplot(q,x,lb,ub,'Normal','cms');

% Print PPCC
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');