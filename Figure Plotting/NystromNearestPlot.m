data = readmatrix('NearestComparionFMM.csv');

filenames = {'NystromFMM.pdf', 'DisjointFMM.pdf', 'IntersectingFMM.pdf'};

for i = 1:3
    plot = reshape(data(i+2,:),5,10);
    figure('Position', [0 0 500 500])
    heatmap(data(1,1:5:end),data(2,1:5),plot,'CellLabelColor','none')
    caxis([0, 1]);
    caxis(log10([min(plot,[],'all') 1]));
    set(gca,'ColorScaling','log') % Log scale
    xlabel("Average force point spacing")
    ylabel('Epsilon')
    set(gca,"FontSize",18)
    savefig(gcf,filenames{i})
end

function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end