classdef Swimmers < handle
    %SWIMMERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0;
        b;
        Geometries;
        swSize;
        swSizeFine;
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
            this.x0 = cell(0);
            this.b = cell(0,2);
            this.Geometries = cell(0);
            this.GeometriesFine = cell(0);
            this.swSize = [];
            this.swSizeFine = [];
            this.NN = cell(0);
            this.Boundary = [];
            this.BoundaryFine = [];
            this.BoundaryNN = [];
        end
        
        function addSwimmer(this,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(varargin) == 3
                this.Geometries{end+1} = varargin{1};
                this.x0{end+1} = varargin{2};
                b = varargin{3};
                this.b(end+1,:) = {b(1:3),b(4:6)};
                this.swSize(end+1) = size(varargin{1},1);
                this.swSizeFine(end+1) = 0;
                this.GeometriesFine{end+1} = [];
                this.NN{end+1} = [];
            elseif length(varargin) == 5
                this.Geometries{end+1} = varargin{1};
                this.swSize(end+1) = size(varargin{1},1);
                this.GeometriesFine{end+1} = varargin{2};
                this.swSizeFine(end+1) = size(varargin{2},1);
                this.NN{end+1} = varargin{3};
                this.x0{end+1} = varargin{4};
                b = varargin{5};
                this.b(end+1,:) = {b(1:3),b(4:6)};
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
            
            noSwimmers = this.swimmerNo();
            [swSizes, swSizesFine] = this.getSwimmerSizes();
            swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
            swIndexFine = zeros(1,noSwimmers+1);swIndexFine(1) = 1;
            for i = 2:noSwimmers+1
                swIndex(i) = swIndex(i-1) + swSizes(i-1);
                swIndexFine(i) = swIndexFine(i-1) + swSizesFine(i-1);
            end
            
            for i = 1:noSwimmers
                b1=this.b{i,1};
                b2=this.b{i,2};
                b3=cross(b1,b2);
                B=[b1(:) b2(:) b3(:)];
                points(swIndex(i):swIndex(i+1)-1,:) = (B * ...
                    points(swIndex(i):swIndex(i+1)-1,:)')';

                points(swIndex(i):swIndex(i+1)-1,:) = ...
                    points(swIndex(i):swIndex(i+1)-1,:)+...
                    this.x0{i}.*ones(swSizes(i),1);
                
                if ~isempty(finePoints)
                    finePoints(swIndexFine(i):swIndexFine(i+1)-1,:) = ...
                        (B * ...
                        finePoints(swIndexFine(i):swIndexFine(i+1)-1,:)')';
                    finePoints(swIndexFine(i):swIndexFine(i+1)-1,:) = ...
                        finePoints(swIndexFine(i):swIndexFine(i+1)-1,:)+...
                        this.x0{i}.*ones(swSizesFine(i),1);
                end
            end
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
            sPoints = vertcat(this.Geometries{:});
            sFinePoints = vertcat(this.GeometriesFine{:});
            sNN = blkdiag(this.NN{:});
            
            bndPoints = this.Boundary;
            bndFinePoints = this.BoundaryFine;
            bndNN = this.BoundaryNN;
            
            points = vertcat(sPoints,bndPoints);
            finePoints = vertcat(sFinePoints,bndFinePoints);
            NN = blkdiag(sNN,bndNN);
            
            noSwimmers = this.swimmerNo();
            [swSizes, swSizesFine] = this.getSwimmerSizes();
            swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
            swIndexFine = zeros(1,noSwimmers+1);swIndexFine(1) = 1;
            for i = 2:noSwimmers+1
                swIndex(i) = swIndex(i-1) + swSizes(i-1);
                swIndexFine(i) = swIndexFine(i-1) + swSizesFine(i-1);
            end
            
            for i = 1:noSwimmers
                b1=this.b{i,1};
                b2=this.b{i,2};
                b3=cross(b1,b2);
                B=[b1(:) b2(:) b3(:)];
                points(swIndex(i):swIndex(i+1)-1,:) = (B * ...
                    points(swIndex(i):swIndex(i+1)-1,:)')';

                points(swIndex(i):swIndex(i+1)-1,:) = ...
                    points(swIndex(i):swIndex(i+1)-1,:)+...
                    this.x0{i}.*ones(swSizes(i),1);
                
                if ~isempty(finePoints)
                    finePoints(swIndexFine(i):swIndexFine(i+1)-1,:) = ...
                        (B * ...
                        finePoints(swIndexFine(i):swIndexFine(i+1)-1,:)')';
                    finePoints(swIndexFine(i):swIndexFine(i+1)-1,:) = ...
                        finePoints(swIndexFine(i):swIndexFine(i+1)-1,:)+...
                        this.x0{i}.*ones(swSizesFine(i),1);
                end
            end
        end
        
        function [swSize,swSizeFine] = getSwimmerSizes(this)
            swSize = this.swSize;
            swSizeFine = this.swSizeFine;
        end
        
        function updateSwimmmer(this,index,x0,b)
            this.x0{index} = x0;
            this.b(index,:) = {b(1:3),b(4:6)};
        end
        
        function this = genTree(this,varargin)
            [points, finePoints, NN] = this.getSwimmerBnd;
            
            this.Tree = OcTree(points,points,'finePoints',finePoints,...
                               'NN',NN,varargin{:});
        end
        
        function this = genKIFMM(this,kernalPar,varargin)
            this.FMM = KIFMM(this.Tree, kernalPar, varargin{:});
        end
    end
end

