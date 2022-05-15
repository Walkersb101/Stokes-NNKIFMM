data = readmatrix('DirectMethodComparison.csv');

figure('Position', [200 200 400 400])
plot(data(1,:),data(3,:)*3,data(1,:),data(6,:)*3)
xlabel('Scalar degrees of freedom')
ylabel('Computation time (seconds)')
legend('KIFMM','Direct Product')
savefig(gcf,'DirectProductCompTime.pdf')

figure('Position', [200 200 400 400])
plot(data(1,:),data(5,:),data(1,:),data(8,:))
set(gca,'yScale','log')
xlabel('Scalar degrees of freedom')
ylabel('Relative error')
legend('KIFMM','Direct Product')
savefig(gcf,'DirectProductComperror.pdf')

function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end