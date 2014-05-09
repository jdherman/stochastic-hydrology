% Function to solve system of MLE equations
function v = gumbel_mle(xtest)
    sebou_data = load('sebou_data.txt');
    x = sebou_data(:,3); % flows, cms
    N = length(x);
    zi = xtest(1);
    alpha = xtest(2);
    v(1) = N - sum(exp(-(x-zi)/alpha));
    v(2) = -N*alpha + sum(x-zi) - sum((x-zi).*exp(-(x-zi)/alpha));
end



