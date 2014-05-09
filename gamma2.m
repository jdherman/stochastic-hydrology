clc; clear all;

sebou_data = load('sebou_data.txt');
x = sebou_data(:,3); % flows, cms
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

ub = gaminv( (i-1)/n + 0.127, ahat,1/bhat);
lb = gaminv( (i)/n - 0.127, ahat,1/bhat);
probplot(q,x,lb,ub,'Gamma-2','cms');