clear all;

data1 = readmatrix('PreconSmall.csv');
data2 = readmatrix('PreconSmallRescale.csv');
data3 = readmatrix('PreconSmallinital.csv');

data3(data3 == 0) = nan;

Normal = data1(1,2:end);
Precon = data1(2,2:end);

ReNormal = data2(1,2:end);
RePrecon = data2(2,2:end);

Inital = data2(1,2:end);

figure('Position', [0 0 1000 1000])
plot(1:numel(Normal),Normal,1:numel(Precon),Precon,...
     1:numel(ReNormal),ReNormal,1:numel(RePrecon),RePrecon,...
     1:numel(Inital),Inital,...
     1:1001,ones(1001,1)*1e-6,'k--');

 
set(gca, 'YScale', 'log')
xlabel('GMRES iteration (Max of 1000)') 
ylabel('Relative Residual error (Target of 1e-6)')
xlim([0,600])

legend('No preconditioner',...
       'Preconditioned',...
       'Rescaled Force, No preconditioner',...
       'Rescaled Force, Preconditioned',...
       'Inital Guess');
   
savefig(gcf,'Rods In Shear Flow Convergence small.pdf')

%-------------------------------------------------------------

data1 = readmatrix('PreconLarge.csv');
data2 = readmatrix('PreconLargeRescale.csv');
data3 = readmatrix('PreconLargeinital.csv');

data3(data3 == 0) = nan;

Normal = data1(1,2:end);
Precon = data1(2,2:end);

ReNormal = data2(1,2:end);
RePrecon = data2(2,2:end);

Inital = data2(1,2:end);

figure('Position', [0 0 1000 1000])
plot(1:numel(Normal),Normal,1:numel(Precon),Precon,...
     1:numel(ReNormal),ReNormal,1:numel(RePrecon),RePrecon,...
     1:numel(Inital),Inital,...
     1:1001,ones(1001,1)*1e-6,'k--');

 
set(gca, 'YScale', 'log')
xlabel('GMRES iteration (Max of 1000)') 
ylabel('Relative Residual error (Target of 1e-6)')
xlim([0,600])

legend('No preconditioner',...
       'Preconditioned',...
       'Rescaled Force, No preconditioner',...
       'Rescaled Force, Preconditioned',...
       'Inital Guess');
   
savefig(gcf,'Rods In Shear Flow Convergence Large.pdf')
   
   
function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end