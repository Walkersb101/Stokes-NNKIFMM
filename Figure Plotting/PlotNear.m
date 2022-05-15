data1 = readmatrix('preconNearDisjoint.csv');
data2 = readmatrix('preconNearJoint.csv');

data2(data2 == 0) = nan;

disNo = data1(1,2:end);
disPre = data1(2,2:end);
joinNo = data2(1,2:end);
joinPre = data2(2,2:end);
figure('Position', [100 100 400 400])
plot(1:numel(disNo),disNo,1:numel(disPre),disPre,...
     1:numel(joinNo),joinNo,1:numel(joinPre),joinPre,...
     1:1001,ones(1001,1)*1e-6,'k--');
set(gca, 'YScale', 'log')
xlabel('GMRES iteration (Max of 1000)') 
ylabel('Relative Residual error (Target of 1e-6)')
xlim([0,1000])

legend('Disjoint',...
'Disjoint Preconditioned',...
'Joint',...
'Joint Preconditioned',...
'Location','northeast','Orientation','vertical');

ylim([1e-7 10])

savefig(gcf,'NearestPrecon.pdf')

function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end
