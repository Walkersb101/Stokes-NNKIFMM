data = readmatrix('PreconSmallEPS.csv');
data(data == 0) = NaN;

figure('Position', [100 100 400 400])
plot(0:1000,data(1,2:end),0:1000,data(2,2:end),0:1000,data(3,2:end),0:1000,ones(1,1001)*1e-6,'--k')
set(gca, 'YScale', 'log')

legend('$\epsilon = 1e-5$: Walltime = 551.613 s',...
       '$\epsilon = 1e-2$: Walltime = 6086.16 s',...
       '$\epsilon = 1e-1$: Walltime = 24050.0 s',...
       'interpreter','latex')

xlim([1 1000])
ylim([1e-7 1])
xticks([21 258 1000])

savefig(gcf,'InitalGuessEPS.pdf')
   
   
function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end