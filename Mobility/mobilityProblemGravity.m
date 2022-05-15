function [out] = mobilityProblemGravity(swimmers,fluidFunc,t,...,
    conditions,kernalPar,TreePar,threads,GPU,blockSize,tol,maxIter,...
    PreconSettings,gravPar)
%mobilityMatrixMulti Solve the matrix inverse for the mobility problem 
%
% Inputs:
%   swimmers       : Swimmers object containing swimmers and boundary
%   fluidFunc      : Function which takes (N,3) array of points and returns
%                    fluid velocity at that point
%   t              : time
%   conditions     : (9*noSwimmers,1) array of (x0;b1;b2)
%   kernalPar      : Kernal parameters
%   TreePar        : Cell array of kernel parameters
%   threads        : Number of threads for KIFMM
%   GPU            : 0 for CPU compute, 1 for GPU compute
%   blockSize      : Size of block matrix computation
%   tol            : Tolerance for GMRES
%   maxIter        : Max number of iterations for GMRES
%   PreconSettings : [Tol,maxIter] for precondioner reduced GMRES. If empty
%                     no preconditioner is used
%
% Outputs:
%   out : (9*noSwimmers,1) of velocity and change in basis vectors

out = zeros(size(conditions));

for i = 0:(numel(conditions)/9)-1
     data = conditions(i*9+1:(i+1)*9);
     swimmers.updateSwimmmer(i+1,data(1:3)',data(4:9))
     swimmers.updateSwimmmerRHS(i+1,gravPar(1:3),gravPar(4));
end

swimmers.genTree(TreePar{:});
swimmers.genKIFMM(kernalPar,'parThreads',threads,'GPU',GPU,...
                  'GMRES',1,'format',2,'blockSize',blockSize);

[B,BT,A1,A2] = genMatrixComponets(swimmers);

pointVel = swimmers.getVelocities();
fluidVel = fluidFunc(t,swimmers.getSwimmers());

vel = pointVel; vel(1:size(fluidVel,1),:) = vel(1:size(fluidVel,1),:) - ...
                                            fluidVel;

[force,moment] = swimmers.getForceMoment();

input = [vel(:);force(1:3:end);force(2:3:end);force(3:3:end);...
                moment(1:3:end);moment(2:3:end);moment(3:3:end)];

if ~isempty(PreconSettings)
    pfunc = @(x) mobilityPrecon(x, swimmers,B,BT,A1,A2,...
                                PreconSettings(1),PreconSettings(2));
    afunc = @(x) mobilityMatrixMulti(pfunc(x),swimmers,B,BT);
            
    x = gmres(afunc, input, [], tol, maxIter);
    x = pfunc(x);
else
    afunc = @(x) mobilityMatrixMulti(x,swimmers,B,BT);
    x = gmres(afunc, input, [], tol, maxIter);
end

x0 = x(end-6*swimmers.swimmerNo()+1:end-3*swimmers.swimmerNo());
omega = x(end-3*swimmers.swimmerNo()+1:end);

x0 = reshape(reshape(x0,[],3).',[],1);
omega = reshape(reshape(omega,[],3).',[],1);

for i = 0:swimmers.swimmerNo()-1
    out((i*9)+1:(i*9)+3) = x0((i*3)+1:(i*3)+3);
    out((i*9)+4:(i*9)+6) = cross(omega((i*3)+1:(i*3)+3),...
        swimmers.b{i+1,1});
    out((i*9)+7:(i*9)+9) = cross(omega((i*3)+1:(i*3)+3),...
        swimmers.b{i+1,2});
end
end