% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

[X,Y,Z] = meshgrid([1e-5,1e-2,1e-1],[1,5,10,20,40], [1,0.7,0.475,0.33]);


for i = 1:numel(X)
    epsilon = X(i);
    noSw = Y(i);
    pointsPerSW = Z(i);

    swimmers = Swimmers();
    swimmersDis = Swimmers();
    swimmersJoin = Swimmers();
    
    for j = 1:noSw
        centre = -20 + 40*rand(1,3);
        points = prolateSpheroid(1,1,pointsPerSW);
        finePoints = prolateSpheroid(1,1,(pointsPerSW/4));
        NN = nearestMatrix(points,finePoints,1,0);
        swimmers.addSwimmer(points,centre,[1;0;0;0;1;0],...
                           zeros(pointsPerSW,3),[0 0 0],[0 0 0])
        swimmersDis.addSwimmer(points,finePoints,NN,centre,...
                           [1;0;0;0;1;0],zeros(pointsPerSW,3),...
                           [0 0 0],[0 0 0])
        finePoints = unique([finePoints;points],'rows');
        NN = nearestMatrix(points,finePoints,1,0);
        swimmersJoin.addSwimmer(points,finePoints,NN,centre,...
                   [1;0;0;0;1;0],zeros(pointsPerSW,3),...
                   [0 0 0],[0 0 0])
    end

    [points,~,~] = swimmers.getSwimmerBnd();
    points = reshape(points.',[],1);

    N = numel(points)/3;

    swimmers.genTree();
    swimmersDis.genTree('nearest',1);
    swimmersJoin.genTree('nearest',1);

    Stokesletmatrix = kernel(points,points,[epsilon,1]);
    [B,BT] = genMatrixComponets(swimmers);

    Mobility = [Stokesletmatrix, full(BT); full(B), zeros(6*noSw)];

    MobilityPre = [eye(3*N),zeros(3*N,6*noSw);...
               full(B)*gather(inv((Stokesletmatrix))), eye(6*noSw)];
           
    B = B./max(B,[],'all'); 
           
    MobilityRescale = [Stokesletmatrix, full(BT); full(B), zeros(6*noSw)];     
    
    [points,fine,NN] = swimmersDis.getSwimmerBnd();
    points = reshape(points.',[],1);
    fine = reshape(fine.',[],1);
    
    StokesletmatrixDis = kernel(fine,points,[epsilon,1])*kron(NN,eye(3));
    [BDis,BTDis] = genMatrixComponets(swimmersDis);

    MobilityDis = [StokesletmatrixDis, full(BTDis); full(BDis), zeros(6*noSw)];

    MobilityDisPre = [eye(3*N),zeros(3*N,6*noSw);...
                    full(BDis)*gather(inv((StokesletmatrixDis))), eye(6*noSw)];
    
    [points,fine,NN] = swimmersJoin.getSwimmerBnd();
    points = reshape(points.',[],1);
    fine = reshape(fine.',[],1);
                
    StokesletmatrixJoin = kernel(fine,points,[epsilon,1])*kron(NN,eye(3));
    [BJoin,BTJoin] = genMatrixComponets(swimmersJoin);

    MobilityJoin = [StokesletmatrixJoin, full(BTJoin); full(BJoin), zeros(6*noSw)];

    MobilityJoinPre = [eye(3*N),zeros(3*N,6*noSw);...
                    full(BJoin)*gather(inv((StokesletmatrixJoin))), eye(6*noSw)];

    filename = sprintf('Condition/%u-%u-%u-%u.csv',...
                        -floor(log10(X(i))),Y(i),...
                        numel(prolateSpheroid(1,1,pointsPerSW))/3,...
                        numel(prolateSpheroid(1,1,pointsPerSW/4))/3);
    
    Stokesletmatrix = calculateData(Stokesletmatrix);
    Mobility = calculateData(Mobility);
    MobilityPre = calculateData(MobilityPre);
    
    MobilityRescale = calculateData(MobilityRescale);
    
    StokesletmatrixDis = calculateData(StokesletmatrixDis);
    MobilityDis = calculateData(MobilityDis);
    MobilityDisPre = calculateData(MobilityDisPre);
    
    StokesletmatrixJoin = calculateData(StokesletmatrixJoin);
    MobilityJoin = calculateData(MobilityJoin);
    MobilityJoinPre = calculateData(MobilityJoinPre);


    padding = numel(MobilityPre)-numel(Stokesletmatrix);
    
    Stokesletmatrix = padarray(Stokesletmatrix,padding,0,'post');
    StokesletmatrixDis = padarray(StokesletmatrixDis,padding,0,'post');
    StokesletmatrixJoin = padarray(StokesletmatrixJoin,padding,0,'post');
    
    file = [Stokesletmatrix,Mobility,MobilityPre,...
            MobilityRescale,...
            StokesletmatrixDis,MobilityDis,MobilityDisPre,...
            StokesletmatrixJoin,MobilityJoin,MobilityJoinPre];
    
    writematrix(file,filename);
                      
    disp(i+"/"+numel(X))
end

function out = calculateData(matrix)
    conditionNum = gather(cond(matrix));
    eigenvalues = gather(eig(matrix));
    matrix = [];
    out = [conditionNum;eigenvalues];
end