Evalres = 50;
[X,Y,Z] = meshgrid(linspace(-3,3,Evalres),linspace(-3,3,Evalres),linspace(-3,3,Evalres));
EvalPoints = [X(:) Y(:) Z(:)];

coarse = prolateSpheroid(0.5,0.5,0.05);
fine = prolateSpheroid(0.5,0.5,0.05/4);
kernelPar = [1e-2,1];


NN = nearestMatrix(coarse,fine,1,0);
vel = reshape(squirmerVel(coarse,2)',[],1);

FMM = KIFMM(OcTree(coarse,coarse,'finePoints',fine,'NN',NN,'nearest',1),...
            kernelPar,'parThreads',4,'GMRES',1);
FMMEval = KIFMM(OcTree(fine,EvalPoints),kernelPar,'parThreads',4);

pot = gmres(@(x) FMM.computeVel(x), vel,[],1e-6,50);

pot = kron(NN,speye(3)) * pot;

evalVel = FMMEval.computeVel(reshape(pot,3,[])');

writematrix(evalVel,"singleFlow.csv")

hold on
[startX,startY,startZ] = meshgrid(-2:0.25:2,-2:0.25:2,2);
verts = stream3(X,Y,Z,reshape(evalVel(:,1),Evalres,Evalres,Evalres),...
                      reshape(evalVel(:,2),Evalres,Evalres,Evalres),...
                      reshape(evalVel(:,3),Evalres,Evalres,Evalres),...
                      startX,startY,startZ);
streamline(verts)
view(3)
plot3(coarse(:,1),coarse(:,2),coarse(:,3),'.');
hold off




function vel = squirmerVel(points,scale)
tangent = [points(:,1).*points(:,3), points(:,2).*points(:,3),...
           -(points(:,1).^2 + points(:,2).^2)]...
           ./(sqrt(points(:,1).^2 + points(:,2).^2+ points(:,3).^2)...
           .*sqrt(points(:,1).^2 + points(:,2).^2));
tangent(isnan(tangent))=0;
theta = acos(points(:,3)./...
        (sqrt(points(:,1).^2 + points(:,2).^2 + points(:,3).^2)));
    
vel = scale*sin(theta).*tangent;
end