clear all;

[X,Y,Z] = meshgrid([1e-5,1e-2,1e-1],[1,5,10,20,40], [1,0.7,0.475,0.33]);

pointarray = arrayfun(@(x) numel(prolateSpheroid(1,1,x)),[1,0.7,0.475,0.33]);

datanames = {{'Stokeslet Matrix'}, {'Mobility Matrix'},...
             {'Mobility Matrix', 'Preconditioned'},...
             {'Mobility Matrix', 'with rescaled total B Matrix'},...
             {'Stokeslet Matrix', 'using Disjoint NEAREST'},...
             {'Mobility Matrix', 'using Disjoint NEAREST'},... 
             {'Mobility Matrix', 'using Disjoint NEAREST', 'Preconditioned'},...
             {'Stokeslet Matrix', 'using Joint NEAREST'},...
             {'Mobility Matrix', 'using Joint NEAREST'},... 
             {'Mobility Matrix', 'using Joint NEAREST', 'Preconditioned'}};

data = cell(size(X));
for i = 1:numel(X)
    filename = sprintf('Condition new/%u-%u-%u-%u.csv',...
                        -floor(log10(X(i))),Y(i),...
                        numel(prolateSpheroid(1,1,Z(i)))/3,...
                        numel(prolateSpheroid(1,1,Z(i)/4))/3);
    data{i} = readmatrix(filename);
end

for j = 1:10
    figure('Position', [100 100 400 400])
    eigenvalues = cell2mat(cellfun(@(x) x(2:end,j),data(5,:,4),'UniformOutput',false));
    for i = 1:3
        hold on;
        plot(sort(nonzeros(abs(eigenvalues(:,i))),'descend'));
        hold off;
    end
    title(datanames{j})
    ylabel({'Magnitude of', 'eigenvalue'})
    xlabel('Eigenvalue Index')
    set(gca, 'YScale', 'log')
    legend('$\epsilon = 10^{-5}$','$\epsilon = 10^{-2}$','$\epsilon = 10^{-1}$','interpreter','latex')
    savefig(gcf,sprintf('Eigen-%s.pdf',strjoin(datanames{j})))
end


function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end