[X,Y,Z] = meshgrid([1e-5,1e-2,1e-1],[1,5,10,20,40], [1,0.7,0.475,0.33]);

pointarray = arrayfun(@(x) numel(prolateSpheroid(1,1,x)),[1,0.7,0.475,0.33]);

datanames = {{'Stokeslet Matrix'}, {'Mobility Matrix'},...
             {'Mobility Matrix', 'Preconditioned'},...
             {'Mobility Matrix', 'with rescaled total B Matrix'},...
             {'Stokeslet Matrix', 'using Disjoint NEAREST'},...
             {'Mobility Matrix', 'using Disjoint NEAREST'},... 
             {'Mobility Matrix', 'using Disjoint NEAREST', 'Preconditioned'},...
             {'Stokeslet Matrix', 'using Contained NEAREST'},...
             {'Mobility Matrix', 'using Contained NEAREST'},... 
             {'Mobility Matrix', 'using Contained NEAREST', 'Preconditioned'}};

data = cell(size(X));
for i = 1:numel(X)
    filename = sprintf('Condition new/%u-%u-%u-%u.csv',...
                        -floor(log10(X(i))),Y(i),...
                        numel(prolateSpheroid(1,1,Z(i)))/3,...
                        numel(prolateSpheroid(1,1,Z(i)/4))/3);
    data{i} = readmatrix(filename);
end

% Res = cell2mat(cellfun(@(x) x(1,1),data(:,:,:),'UniformOutput',false));
% Mob = cell2mat(cellfun(@(x) x(1,2),data(:,:,:),'UniformOutput',false));
% MobPre = cell2mat(cellfun(@(x) x(1,3),data(:,:,:),'UniformOutput',false));
% MobN = cell2mat(cellfun(@(x) x(1,4),data(:,:,:),'UniformOutput',false));
% MobPreN = cell2mat(cellfun(@(x) x(1,5),data(:,:,:),'UniformOutput',false));
% MobR = cell2mat(cellfun(@(x) x(1,6),data(:,:,:),'UniformOutput',false));
% MobPreR = cell2mat(cellfun(@(x) x(1,7),data(:,:,:),'UniformOutput',false));
% MobNR = cell2mat(cellfun(@(x) x(1,8),data(:,:,:),'UniformOutput',false));
% MobPreNR = cell2mat(cellfun(@(x) x(1,9),data(:,:,:),'UniformOutput',false));

figure('Position', [200 200 600 600])
% t = tiledlayout(5,2,'TileSpacing','compact');
for i = 1:10
%     ax = nexttile;
    plot = cell2mat(cellfun(@(x) x(1,i),data(:,:,:),'UniformOutput',false));
    plot = permute(plot(:,1,:),[1 3 2]);
    heatmap(pointarray(1:4),[1,5,10,20,40],plot,'CellLabelColor','none','FontSize',16,'PositionConstraint','outerposition');
%     title(datanames{i})
    xlabel({'Scalar degrees of', 'freedom per swimmer'})
    ylabel('Number of Swimmers')
    savefig(gcf,sprintf('%s-%u.pdf',strjoin(datanames{i}),5))
end

figure('Position', [200 200 600 600])
% t = tiledlayout(5,2,'TileSpacing','compact');
for i = 1:10
%     ax = nexttile;
    plot = cell2mat(cellfun(@(x) x(1,i),data(:,:,:),'UniformOutput',false));
    plot = permute(plot(:,2,:),[1 3 2]);
    heatmap(pointarray(1:4),[1,5,10,20,40],plot,'CellLabelColor','none','FontSize',16);
%     title(datanames{i})
    xlabel({'Scalar degrees of', 'freedom per swimmer'})
    ylabel('Number of Swimmers')
    savefig(gcf,sprintf('%s-%u.pdf',strjoin(datanames{i}),2))
end

figure('Position', [200 200 600 600])
% t = tiledlayout(5,2,'TileSpacing','compact');
for i = 1:10
%     ax = nexttile;
    plot = cell2mat(cellfun(@(x) x(1,i),data(:,:,:),'UniformOutput',false));
    plot = permute(plot(:,3,:),[1 3 2]);
    heatmap(pointarray(1:4),[1,5,10,20,40],plot,'CellLabelColor','none','FontSize',16);
%     title(datanames{i})
    xlabel({'Scalar degrees of', 'freedom per swimmer'})
    ylabel('Number of Swimmers')
    savefig(gcf,sprintf('%s-%u.pdf',strjoin(datanames{i}),1))
end



function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)*1.1, pos(4)*1.1])
print(fig,name,'-dpdf','-r0')
end