function ax = newsubplot(position, xlab, ylab)
    % Creates new subplot in specified position on current figure
    % with xlab xlabel and ylab ylabel
    ax = subplot(position); 
    hold on
    set(ax,'FontSize',7) %and other properties
    xlabel(xlab);
    ylabel(ylab,'interpreter','latex');
    grid on;
end