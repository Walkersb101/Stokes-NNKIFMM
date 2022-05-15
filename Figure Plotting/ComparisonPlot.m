CPU = readmatrix('DirectMethodComparison');
GPU = readmatrix('DirectMethodComparisonGPU');


figure('Position', [200 200 400 400])

hold on;

plot(CPU(1,:),CPU(3,:)*3,CPU(1,:),CPU(6,:)*3,...
     GPU(1,:),GPU(3,:)*3,GPU(1,:),GPU(6,:)*3);
 
plot(CPU(1,:),3e-6*(CPU(1,:)/3).^2,'k--', CPU(1,:),1.66e-2*(CPU(1,:)- 5e4),'k-.')
xlabel('Scalar degrees of freedom')
ylabel('Computation time (seconds)')
legend('KIFMM CPU','Direct Product CPU', 'KIFMM GPU','Direct Product GPU','$\mathcal{O}(N^2)$','$\mathcal{O}(N)$','interpreter','latex','Location','northwest')
% ylim(log10([0 max(CPU(6,:)*3)*1.1]))
box on;
savefig(gcf,'DirectProductCompTime.pdf')

figure('Position', [200 200 400 400])
plot(CPU(1,:),CPU(5,:),CPU(1,:),CPU(8,:),...
     GPU(1,:),GPU(5,:),GPU(1,:),GPU(8,:));
set(gca,'yScale','log')
xlabel('Scalar degrees of freedom')
ylabel('Relative error')
legend('KIFMM CPU','Direct Product CPU', 'KIFMM GPU','Direct Product GPU','interpreter','latex')
box on;
savefig(gcf,'DirectProductComperror.pdf')

function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end