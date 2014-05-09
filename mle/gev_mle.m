% Function to solve system of MLE equations
function v = gev_mle(xtest)
    sebou_data = load('sebou_data.txt');
    x = sebou_data(:,3); % flows, cms
    
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
    v = -1*(-N*log(a) + sum((1/k-1)*log(y) - y.^(1/k)));
end



