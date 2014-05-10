function [pe,pet] = censoredpp(x,xtin)
    
    % Assign plotting probabilities to censored data
    % Inputs:
    % Observed data x (e.g. [14,21,35,43,52])
    % Censored data xtin (e.g. [10 10 10 20 30]) - Repeat values needed!
    % Outputs:
    % pe: estimated exceedance probabilities for observed data
    % (column 1: data value, column 2: pe)
    % pet: estimated exceedance probs. for thresholds (not used for plot?)
    
    % NOTE: THIS ASSUMES CENSORED DATA ARE "LESS THAN" VALUES
    
    N = length(xtin) + length(x);
    xt = unique(xtin);
    
    xt = sort(xt, 'descend');
    x = sort(x, 'descend');
    pet = zeros(length(xt),1); % exceedance prob. for thresholds
    pe = zeros(length(x),1); % exceedance prob. for observations
    
    % First, assign exceedance prob. to thresholds
    % as they are assigned, subtract count off of N
    % The first one is a little different than the others
    pet(1) = length(x(x > xt(1)))/N;
    N = N - length(x(x > xt(1))) - length(xtin(xtin == xt(1)));
    
    for i=2:length(xt)
        pet(i) = pet(i-1) + (1-pet(i-1))*length(x(x > xt(i) & x < xt(i-1)))/N;
        N = N - length(x(x > xt(i) & x < xt(i-1))) - length(xtin(xtin == xt(i)));
    end
    
    % Then assign exceedance prob. to observations
    % Need to grab those between each threshold
    xidx = find(x > xt(1));
    for j=1:length(xidx)
        pe(xidx(j)) = j/(length(xidx)+1)*pet(1);
    end
    
    for i=2:length(xt)
        xidx = find(x < xt(i-1) & x > xt(i));

        for j=1:length(xidx)
            pe(xidx(j)) = pet(i-1) + j/(length(xidx)+1)*(pet(i) - pet(i-1));
        end
    end
    
    pe = [x',pe];
    pet = [xt',pet];
    
end