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
    for i=1:length(xt)
        if(i==1)
            xs = x(x > xt(i));
            pet(i) = length(xs)/N;
        else
            xs = x(x > xt(i) & x < xt(i-1));
            pet(i) = pet(i-1) + (1-pet(i-1))*length(xs)/N;
        end
        N = N - length(xs) - length(xtin(xtin == xt(i)));
        
        % Then assign probs. to observations that fall between thresholds
        for j=1:length(xs)
            idx = find(x==xs(j));
            if(i==1)
                pe(idx) = j/(length(xs)+1)*pet(i);
            else
                pe(idx) = pet(i-1) + j/(length(xs)+1)*(pet(i) - pet(i-1));
            end
        end
    end
    
    % return as nonexceedance for probability plots
    pe = flipud([x', 1-pe]);
    pet = flipud([xt', 1-pet]);
    
end