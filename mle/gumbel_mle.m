% Function to solve system of MLE equations
function v = gumbel_mle(xtest)
    addpath ..
    x = load('input_data.txt');
    N = length(x);
    zi = xtest(1);
    alpha = xtest(2);
    v(1) = N - sum(exp(-(x-zi)/alpha));
    v(2) = -N*alpha + sum(x-zi) - sum((x-zi).*exp(-(x-zi)/alpha));
end



