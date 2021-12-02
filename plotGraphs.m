file = load('ComputeTimeGPU2');
storage = file.temp;



ax1 = subplot(2,1,1);
x = linspace(min(storage(1,4:end)),max(storage(1,:)),50);
plot(ax1,storage(1,3:end),storage(2,3:end),'b',storage(1,3:end),storage(3,3:end),'r',x,(1.5e-9)*x.^2,'k--',x,(1.4e-5)*x.*log(x),'k-.');

ax1.XScale = 'log';
ax1.YScale = 'log';

xlabel('Degrees of freedom') 
ylabel('Computational time [s]') 

xlim([min(storage(1,3:end))*0.9 max(storage(1,3:end))*1.1])

legend('Direct solver','KIFMM', 'N^2', 'NlogN','Location','northwest')

ax2 = subplot(2,1,2);
plot(storage(1,3:end),storage(4,3:end),'k')

ax2.XScale = 'log';
ax2.YScale = 'log';

xlabel('Degrees of freedom') 
ylabel('Relative Error between direct and KIFMM solution') 

xlim([min(storage(1,3:end))*0.9 max(storage(1,3:end))*1.1])