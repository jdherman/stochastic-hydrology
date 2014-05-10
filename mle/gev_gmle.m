% Function to solve system of MLE equations
function v = gev_gmle(xtest)
    addpath ..
    x = load('input_data.txt');
    
    N = length(x);
    
    a = xtest(1);
    k = xtest(2);
    q = xtest(3);

    % enforce constraints on distribution
    % search over q but evaluate with z
    if(k < 0)
        z = min(x) - a/k - exp(q);
    else
        z = max(x) - a/k + exp(q);
    end
    
    % transformed y
    y = 1 - (k/a)*(x-z); 
    
    p = 6;
    q = 9;
    v = -1*(-N*log(a) + sum((1/k-1)*log(y) - y.^(1/k)) + ...
        (p-1)*log(0.5+k) + (q-1)*log(0.5-k) - log(beta(p,q)));
    
end



