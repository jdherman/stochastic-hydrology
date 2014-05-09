clc; clear all;

sebou_data = load('sebou_data.txt');
x = sebou_data(:,3); % flows, cms

% Normal distribution
m = mean(x); v = var(x);
disp('Normal distribution (MM/MLE)');
disp(['m = ' num2str(m) ', s2 = ' num2str(v) ' (s = ' num2str(sqrt(v)) ')']);
disp(['1%: ' num2str(norminv(.01,m,sqrt(v))) ', 99%: ' num2str(norminv(.99,m,sqrt(v))) ' m^3/s']);

% Probability Plot w/90% KS Bounds
i = 1:length(x);
n = length(x);
q = norminv( (i-3/8)/(n+1/4) , m, sqrt(v));
ub = norminv( (i-1)/n + 0.127, m, sqrt(v));
lb = norminv( (i)/n - 0.127, m, sqrt(v));
probplot(q,x,lb,ub,'Normal','cms');

% Print PPCC
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');