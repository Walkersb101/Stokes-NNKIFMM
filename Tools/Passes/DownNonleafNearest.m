function pot = DownNonleafNearest(vsurf,vpot,xsurf,xpot,parentsurf,parentpot,upsurf,downsurf,NN,blockSize,kernelPar,GPU)
% UpNonleafNearest  Compute the affect of nodes on the downward equivlient 
% surface for non leaf nodes
% Compute the efffect of points on the upward equivilent potential
% before inverting to find the equilent potential on the downward surface
%Inputs:
%   vsurf      : A vector of positions at v nodes upward surfaces
%   vpot       : A vector of potentials at v nodes upward surfaces
%   xsurf      : A vector of positions at v nodes downward surfaces
%   xpot       : A vector of potentials at v nodes downward surfaces
%   parentsurf : A vector of positions at the parent node downward surface
%   parentpot  : A vector of potentials at the parent node downward 
%                surfaces
%   upsurf     : A (3*Nq,1) vector defining the position of the upward 
%                surface
%   downsurf   : A (3*Nq,1) vector defining the position of the 
%                downward surface
%   NN         : (Q, N) Matrix for interpolation, empty for non nearest
%                application
%   kernelPar  : Paramaters for the kernel
%   GPU        : 1 for GPU compute, 0 otherwise
%
%Output:
%   pot : (1, 3*Nq) vector of potentials at the upward surface

if GPU
    vsurf = gpuArray(vsurf);
    vpot = gpuArray(vpot);
    xsurf = gpuArray(xsurf);
    xpot = gpuArray(xpot);
    parentsurf = gpuArray(parentsurf);
    parentpot = gpuArray(parentpot);
    upsurf = gpuArray(upsurf);
    downsurf = gpuArray(downsurf);
    NN = gpuArray(NN);
    
    RHS = blockcomputation(vsurf,upsurf,vpot,blockSize,kernelPar) + ...
      nearestBlockComputation(xsurf,upsurf,NN,xpot,blockSize,kernelPar) + ...
      blockcomputation(parentsurf,upsurf,parentpot,blockSize,kernelPar);

    pot = gather((kernel(downsurf,upsurf,kernelPar) \ RHS)');
    
    clear vsurf vpot xsurf xpot parentsurf parentpot upsurf downsurf RHS;
else
    RHS = blockcomputation(vsurf,upsurf,vpot,blockSize,kernelPar) + ...
      nearestBlockComputation(xsurf,upsurf,NN,xpot,blockSize,kernelPar) + ...
      blockcomputation(parentsurf,upsurf,parentpot,blockSize,kernelPar);
  
    pot = (kernel(downsurf,upsurf,kernelPar) \ RHS)';
end
end

