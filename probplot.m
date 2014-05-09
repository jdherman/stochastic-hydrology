function ax = probplot(q,x,lb,ub,distname,units)
    plot(q,q, 'linewidth', 2, 'color', 0.75*ones(1,3));
    hold on;
    plot(q,ub, '--', 'linewidth', 2, 'color', 0.75*ones(1,3));
    plot(q,lb, '--', 'linewidth', 2, 'color', 0.75*ones(1,3));
    plot(q, sort(x), '.k', 'markersize', 13);
    title([distname ' Q-Q Plot with 90% KS Bounds']);
    xlabel(['Quantiles of Fitted ' distname ' Distribution (' units ')']);
    ylabel(['Observed Values X(i) and KS Bounds (' units ')']);
    ax = gca;
end