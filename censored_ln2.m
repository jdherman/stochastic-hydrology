clc; clear all;
addpath mle;

% 2-parameter Lognormal distribution for censored data
x = [55 72 186 94 64 142 48 37 75 29 31 54]; % observations
xt = [50 50 50 50 50 50 25 25 25]; % "below threshold"
[pe,pet] = censoredpp(x,xt);

% Probability plot regression
disp('LN2 Prob. Plot Regression');
plot(logninv(pe(:,2)),log(sort(x)),'o');
title('Probability plot of censored data');
xlabel('Quantiles of standard lognormal distribution');
ylabel('Log(x) observations');
p = polyfit(logninv(pe(:,2)),log(sort(x')),1);
hold on; lsline;
disp(['m = ' num2str(p(2)) ', s2 = ' num2str(p(1)^2) ' (s = ' num2str(p(1)) ')']);
disp(' ');

% MLE for censored data
% This is a numerical problem because of the censored points
disp('LN2 Censored MLE');
disp('WARNING: double-check x data in the MLE function');
optpt = fminsearch(@ln2_censored_mle, [p(1) p(2)]);
m = optpt(1);
s = optpt(2);
disp(['m = ' num2str(m) ', s2 = ' num2str(s^2) ' (s = ' num2str(s) ')']);
