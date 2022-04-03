classdef Swimmers < handle
    %SWIMMERS Struture to hold data on swimmers and boundary for mobility
    % problems. Swimmer geometry is stored in Cell arrays, If Nearest
    % neibour interpotation is selected the fine point array and nearest
    % neigbour matrix is stored in seperated cell arrays.
    %
    % Properties:
    %   x0             : center of each swimmer
    %   b              : Basic of each swimmer, first two stored with 3rd
    %                    calculated from their cross product
    %   Geometries     : Storage of force points for each swimmer
    %   Vel            : Storage of velocity atforce points for each 
    %                    swimmer
    %   swSize         : Number of force points for each swimmer
    %   GeometriesFine : Storage of fine quadrature points
    %   swSizeFine     : Number of fine quadrature points for each swimmer
    %   NN             : Nearest neigbour matrix for each swimmer
    %   F              : Total force acting on each swimmer 
    %   Moment         : Total moment acting on each swimmer
    %   Boundary       : Location of force points on the boundary
    %   BoundaryVel    : Velocity of boundary at each force point
    %   BoundaryFine   : Location of fine qudrature rule on the boudnary
    %   BoundaryNN     : Nearest neigbour matrix for the boundary
    %   Tree           : Tree for KIFMM computation
    %   FMM            : FMM object for KIFMM computation
    %
    % Functions:
    %   addSwimmer      : Add Swimmer to problem
    %   updateBnd       : Add boundary to problem
    %   getSwimmers     : Get point data from all swimmers
    %   getIndSwimmer   : Get data from individual swimmer
    %   swimmerNo       : Returns number of swimmers
    %   getBnd          : Get data about Boundary
    %   getSwimmersBnd  : Get data from swimmers and boundary
    %   getVelocities   : Extract velocities from swimmers and boundaries
    %   getSwimmerSizes : Return number of points for each swimmer
    %   getForceMoment  : Extract total force and moment on each body
    %   updateSwimmmer  : Update the centre and basis of each swimmer
    %   genTree         : Generate adaptive Octree for the computation the 
    %                     KIFMM method
    %   genKIFMM        : Generate KIFMM handler for the computation the 
    %                     KIFMM method
    
    properties
        x0;
        b;
        Geometries;
        Vel;
        swSize;
        GeometriesFine;
        swSizeFine;
        NN;
        F;
        Moment;
        Boundary;
        BoundaryVel;
        BoundaryFine;
        BoundaryNN;
        Tree;
        FMM;
    end
    
    methods
        function this = Swimmers()
            %SWIMMERS Initialsie empty arrays to store bodies for mobility
            % problems
            this.x0 = cell(0);
            this.b = cell(0,2);
            this.Geometries = cell(0);
            this.Vel = cell(0);
            this.GeometriesFine = cell(0);
            this.swSize = [];
            this.swSizeFine = [];
            this.NN = cell(0);
            this.F = cell(0);
            this.Moment = cell(0);
            this.Boundary = [];
            this.BoundaryVel = [];
            this.BoundaryFine = [];
            this.BoundaryNN = [];
        end
        
        function addSwimmer(this,varargin)
            %addSwimmer add Swimmer to problem
            %   Adds swimmer to arrays, if standard problem then the number
            %   of inputs is 6 structured
            %   (points,centre,basis,vel,F,Moment) 
            %   If NEAREST is choosen then 8 inputs are required as 
            %   (points,fine points,NN,centre,basis,vel,F,Moment).
            %
            %   Points and finePoints are structed as a (N,3) and (Q,3)
            %   array respectivly, NN is a (Q,N) matrix.
            %   center is a (1,3) array, basis is a (6,1) array of the
            %   first two basis vectors (b1;b2), vel is a (N,3) array of
            %   points. F and Moment are (3,1)vectors.
            
            if length(varargin) == 6
                this.Geometries{end+1} = varargin{1};
                this.x0{end+1} = varargin{2};
                b = varargin{3};
                this.b(end+1,:) = {b(1:3),b(4:6)};
                this.swSize(end+1) = size(varargin{1},1);
                this.swSizeFine(end+1) = 0;
                this.GeometriesFine{end+1} = [];
                this.NN{end+1} = [];
                this.Vel{end+1} = varargin{4};
                this.F{end+1} = varargin{5};
                this.Moment{end+1} = varargin{6};
            elseif length(varargin) == 8
                this.Geometries{end+1} = varargin{1};
                this.swSize(end+1) = size(varargin{1},1);
                this.GeometriesFine{end+1} = varargin{2};
                this.swSizeFine(end+1) = size(varargin{2},1);
                this.NN{end+1} = varargin{3};
                this.x0{end+1} = varargin{4};
                b = varargin{5};
                this.b(end+1,:) = {b(1:3),b(4:6)};
                this.Vel{end+1} = varargin{6};
                this.F{end+1} = varargin{7};
                this.Moment{end+1} = varargin{8};
            else
                warning("Invalid number of inputs")
            end

        end
        
        function updateBnd(this,varargin)
            %updateBnd Create boundary for problem
            %   Create boundary for problem, if number of inputs is 2 then
            %   normal geomertry and velocity are stored as (N,3) arrays.
            %   NEAREST method requries 4 inputs 
            %   (points,fine points, NN, vel). 
            %   Fine points is a (Q,3) array and NN is a (Q,N) matrix.
            
            if length(varargin) == 2
                this.Boundary = varargin{1};
                this.BoundaryFine = double.empty(0,3);
                this.BoundaryNN = double.empty(0,0);
                this.BoundaryVel = varargin{2};
            elseif length(varargin) == 4
                this.Boundary = varargin{1};
                this.BoundaryFine = varargin{2};
                this.BoundaryNN = varargin{3};
                this.BoundaryVel = varargin{4};
            else
                warning("Invalid number of inputs")
            end
        end
        
        function [points, finePoints, NN, pointsrot] = getSwimmers(this)
            %getSwimmers Get data from swimmers
            %   Join data from all swimmers into a single arrays which can
            %   be used for computation. 
            %
            % Output:
            %   points     : rotated and translated force points (N,3)
            %   finePoints : rotated and translated fine quadrature points
            %                (Q,3). Empty array is returned if non-NEAREST
            %   NN         : Joined NEAREST Neighbour matrix. 
            %                Empty array is returned if non-NEAREST
            %   pointsrot  : Rotated force points. NON TRANSLATED.
            
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
                pointsrot(swIndex(i):swIndex(i+1)-1,:) = (B * ...
                    points(swIndex(i):swIndex(i+1)-1,:)')';

                points(swIndex(i):swIndex(i+1)-1,:) = ...
                    pointsrot(swIndex(i):swIndex(i+1)-1,:)+...
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
        
        function [points, finePoints, NN, x0] = getIndSwimmer(this,index)
            %getIndSwimmer Get data from individual swimmer
            %   extract all data for an individual swimmer
            %
            % Input:
            %   index : interger index of swimmer
            %
            % Output:
            %   points     : rotated and translated force points (N,3)
            %   finePoints : rotated and translated fine quadrature points
            %                (Q,3). Empty array is returned if non-NEAREST
            %   NN         : Joined NEAREST Neighbour matrix. 
            %                Empty array is returned if non-NEAREST
            %   x0         : center of swimmer 
            
            points = this.Geometries{index};
            finePoints = this.GeometriesFine{index};
            NN = this.NN{index};
            
            b1=this.b{index,1};
            b2=this.b{index,2};
            b3=cross(b1,b2);
            B=[b1(:) b2(:) b3(:)];
            points = (B * points')';
            if ~isempty(finePoints)
                finePoints = (B * finePoints')';
            end
            x0 = this.x0{index};
        end
        
        function no = swimmerNo(this)
            %swimmerNo returns number of swimmers
            
            no = length(this.Geometries);
        end
        
        function [points, finePoints, NN] = getBnd(this)
            %getBnd Get data about Boundary
            %   Extract data about boundary
            %
            % Output:
            %   points     : force points 
            %   finePoints : fine quadrature points
            %                Empty array is returned if non-NEAREST
            %   NN         : NEAREST Neighbour matrix. 
            %                Empty array is returned if non-NEAREST
            
            points = this.Boundary;
            finePoints = this.BoundaryFine;
            NN = this.BoundaryNN;
        end
        
        function [points, finePoints, NN] = getSwimmerBnd(this)
            %getSwimmersBnd Get data from swimmers and boundary
            %   Join data from all swimmers and boundary into a 
            %   single arrays which can be used for computation. 
            %   Swimmers are placed first followed by boundary.
            %
            % Output:
            %   points     : Force points (N,3)
            %   finePoints : Fine quadrature points (Q,3)
            %                Empty array is returned if non-NEAREST
            %   NN         : Joined NEAREST Neighbour matrix. 
            %                Empty array is returned if non-NEAREST

            
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
        
        function Vel = getVelocities(this)
            %getVelocities Extract velocities from swimmers and boundaries
            %   rotate velocites from body frame for swimmers and append
            %   velocities of boundary 
            %
            % Output:
            %   Vel     : velocities at force points
            
            swVel = vertcat(this.Vel{:});
            bndVel = this.BoundaryVel;
            
            noSwimmers = this.swimmerNo();
            [swSizes, ~] = this.getSwimmerSizes();
            swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
            for i = 2:noSwimmers+1
                swIndex(i) = swIndex(i-1) + swSizes(i-1);
            end
            
            for i = 1:noSwimmers
                b1=this.b{i,1};
                b2=this.b{i,2};
                b3=cross(b1,b2);
                B=[b1(:) b2(:) b3(:)];
                swVel(swIndex(i):swIndex(i+1)-1,:) = (B * ...
                swVel(swIndex(i):swIndex(i+1)-1,:)')';
            end
            
            Vel = vertcat(swVel,bndVel);
        end
        
        function [Force,Moment] = getForceMoment(this)
            %getForceMoment Extract total force and moment on each body
            %   calcualte force and moment on each body in labratory frame
            %
            % Output:
            %   Force  : Force on each swimmers
            %   Moment : Moment on each swimmer
            
            Force = vertcat(this.F{:});
            Moment = vertcat(this.Moment{:});
            
            noSwimmers = this.swimmerNo();
            for i = 0:noSwimmers-1
                b1=this.b{i+1,1};
                b2=this.b{i+1,2};
                b3=cross(b1,b2);
                B=[b1(:) b2(:) b3(:)];
                Force((3*i)+1:(3*i)+3) = (B * Force((3*i)+1:(3*i)+3))';
                Moment((3*i)+1:(3*i)+3) = (B * Moment((3*i)+1:(3*i)+3))';
            end
            
        end
        
        function [swSize,swSizeFine] = getSwimmerSizes(this)
            %getSwimmerSizes return number of points for each swimmer
            
            swSize = this.swSize;
            swSizeFine = this.swSizeFine;
        end
        
        function updateSwimmmer(this,index,x0,b)
            %updateSwimmmer update the centre and basis of each swimmer
            %
            % Input :
            %   x0  : (1,3) array denoting the centre of the array
            %   b   : (6,1) array denoting the first two basis vectors
            
            this.x0{index} = x0;
            this.b(index,:) = {b(1:3),b(4:6)};
        end
        
        function this = genTree(this,varargin)
            %genTree Generate adaptive Octree for the computation the KIFMM
            % method
            %
            % Input : all inputs passed to Octree class, see Octree class
            %         for details
            
            [points, finePoints, NN] = this.getSwimmerBnd;
            
            this.Tree = OcTree(points,points,'finePoints',finePoints,...
                               'NN',NN,varargin{:});
        end
        
        function this = genKIFMM(this,kernalPar,varargin)
            %genKIFMM Generate KIFMM handler for the computation the KIFMM
            % method
            %
            % Input : 
            %   kernalPar : parapeters passed to KIFMM kernel
            %   other     : all other inputs passed to Octree class, 
            %         see Octree classfor details
            
            this.FMM = KIFMM(this.Tree, kernalPar, varargin{:});
        end
    end
end
