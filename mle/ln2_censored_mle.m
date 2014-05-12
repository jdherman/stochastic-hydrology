% Function to solve system of MLE equations
function v = ln2_censored_mle(xtest)
    
    x = [55 72 186 94 64 142 48 37 75 29 31 54]; % observations
    xt = [50 50 50 50 50 50 25 25 25]; % "below threshold"
    m = xtest(1);
    s = xtest(2);
    
    L = 1;
    for i=1:length(xt)
        L = L*logncdf(xt(i),m,s);
    end
    
    for i=1:length(x)
        L = L*lognpdf(x(i),m,s);
    end
    
    v(1) = -L;
end



