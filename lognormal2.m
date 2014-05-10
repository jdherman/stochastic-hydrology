clc; clear all;

x = load('input_data.txt');
i = 1:length(x);
n = length(x);

% 2-parameter lognormal: MM (LB Eqn 7.66)
disp('2-parameter lognormal distribution, MM');
vhat = log(1 + var(x)/(mean(x)^2));
mhat = log(mean(x)) - (1/2)*vhat;
disp(['m = ' num2str(mhat) ', s2 = ' num2str(vhat) ' (s = ' num2str(sqrt(vhat)) ')']);
disp(['1%: ' num2str(logninv(.01,mhat,sqrt(vhat))) ', 99%: ' num2str(logninv(.99,mhat,sqrt(vhat))) ' m^3/s']);
q = logninv( (i-3/8)/(n+1/4) , mhat, sqrt(vhat));
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% 2-parameter lognormal: MLE (LB 7.63)
disp('2-parameter lognormal distribution, MLE');
mhat = mean(log(x));
vhat = mean((log(x) - mhat).^2);
disp(['m = ' num2str(mhat) ', s2 = ' num2str(vhat) ' (s = ' num2str(sqrt(vhat)) ')']);
disp(['1%: ' num2str(logninv(.01,mhat,sqrt(vhat))) ', 99%: ' num2str(logninv(.99,mhat,sqrt(vhat))) ' m^3/s']);
q = logninv( (i-3/8)/(n+1/4) , mhat, sqrt(vhat));
r = corrcoef(sort(x),q);
disp(['PPCC: ' num2str(r(1,2))]);
disp(' ');

% Probability plot for MLE (or MM if you cut and paste)
ub = logninv( (i-1)/n + 0.127, mhat, sqrt(vhat));
lb = logninv( (i)/n - 0.127, mhat, sqrt(vhat));
probplot(q,x,lb,ub,'Lognormal-2','cms');

