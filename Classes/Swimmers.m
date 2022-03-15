classdef Swimmers < handle
    %SWIMMERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0;
        Geometries;
        swSize;
        GeometriesFine
        NN;
        Boundary;
        BoundaryFine;
        BoundaryNN;
        Tree;
        FMM;
    end
    
    methods
        function this = Swimmers()
            %SWIMMERS Construct an instance of this class
            %   Detailed explanation goes here
            this.Geometries = cell(0);
            this.GeometriesFine = cell(0);
            this.swSize = [];
            this.NN = cell(0);
            this.Boundary = [];
            this.BoundaryFine = [];
            this.BoundaryNN = [];
        end
        
        function addSwimmer(this,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(varargin) == 2
                this.Geometries{end+1} = varargin{1};
                this.x0{end+1} = varargin{2};
                this.swSize(end+1) = size(varargin{1},1);
                this.GeometriesFine{end+1} = [];
                this.NN{end+1} = [];
            elseif length(varargin) == 4
                this.Geometries = varargin{1};
                this.swSize(end+1) = size(varargin{1},1);
                this.GeometriesFine{end+1} = varargin{2};
                this.NN{end+1} = varargin{3};
                this.x0{end+1} = varargin{4};
            else
                warning("Invalid number of inputs")
            end

        end
        
        function updateBnd(this,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(varargin) == 1
                this.Boundary = varargin{1};
                this.BoundaryFine = 0;
                this.BoundaryNN =0;
            elseif length(varargin) == 3
                this.Boundary = varargin{1};
                this.BoundaryFine = varargin{2};
                this.BoundaryNN = varargin{3};
            else
                warning("Invalid number of inputs")
            end
        end
        
        function [points, finePoints, NN] = getSwimmers(this)
            points = vertcat(this.Geometries{:});
            finePoints = vertcat(this.GeometriesFine{:});
            NN = blkdiag(this.NN{:});
        end
        
        function [points, finePoints, NN] = getIndSwimmer(this,index)
            points = this.Geometries{index};
            finePoints = this.GeometriesFine{index};
            NN = this.NN{index};
        end
        
        function no = swimmerNo(this)
            no = length(this.Geometries);
        end
        
        function [points, finePoints, NN] = getBnd(this)
            points = this.Boundary;
            finePoints = this.BoundaryFine;
            NN = this.BoundaryNN;
        end
        
        function [points, finePoints, NN] = getSwimmerBnd(this)
            [sPoints,sFinePoints,sNN] = this.getSwimmers();
            [bndPoints,bndFinePoints,bndNN] = this.getBnd();
            points = vertcat(sPoints,bndPoints);
            finePoints = vertcat(sFinePoints,bndFinePoints);
            NN = blkdiag(sNN,bndNN);
        end
        
        function swSize = getSwimmerSizes(this)
            swSize = this.swSize;
        end
        
        function this = genTree(this,varargin)
            [points, finePoints, NN] = this.getSwimmerBnd;
            this.Tree = OcTree(points,points,'finePoints',finePoints,'NN',NN,varargin{:});
        end
        
        function this = genKIFMM(this,kernalPar,varargin)
            this.FMM = KIFMM(this.Tree, kernalPar, varargin{:});
        end
    end
end

