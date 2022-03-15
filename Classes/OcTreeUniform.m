classdef OcTreeUniform < handle
    % OCTree Compute OcTree for input points
    %   Generate OcTree for input points. The root node will be the
    %   smallest cube which encompases all points given. The this will
    %   adaptivly divide until either the max depth is reaches or all
    %   leaves have less than the maximum capacity. Data is stored
    %   linearly with each row representing a new node. 
    %
    % Inputs:
    %   potPoints    : (N,3) array storing the postion of each potential 
    %                   point
    %   targetPoints : (M,3) array of target points
    %   arguement    : See Optional Arguments below
    %
    % Optional Arguments:
    %   nodeCapacity : sets the maximum capacity of the node(default = 200)
    %   maxDepth     : sets the maximum number of node divsions
    %                  (default = 21)
    %   nearest      : 0 for normal, 1 for Nearest interpolation
    %   finepoints   : (Q,3) array of fine quadrature points for nearest
    %   NN           : (Q, N) sparce Matrix for interpolation, generated if
    %                   not given

    %
    % Properties:
    %   points       : (N,3) array storing the postion of each point
    %   potentials   : (N,1) logic array to slice points into potential
    %                  points
    %   targets      : (N,1) logic array to slice points into target
    %                  points
    %   finePoints   : (Q,3) array of finerature points)
    %   NN           : (Q,N) sparce nearest interpolation matrix
    %   pointIndex   : Store the node index in which the point is in
    %   nodeCorners  : Stores the coorinates of the corners of the node
    %                  [x1 y1 z1 x2 y2 z2] where point 1 is the lower left
    %                  front corner and point 2 is the upper right back
    %                  corner.
    %   nodeLevel    : Stores which level the node is on 
    %                  (0 is the root node)
    %   nodeParent   : Stores the index of the parent of the current node
    %   nodeChildren : Stores the index of the children of the node as 
    %                  (1,8) array
    %   nodeCount    : Stores the total number of nodes
    %   nodeCapacity : Stores the maximum number of points in a bin
    %
    % Functions:
    %   OcTree           : Initialsied the OcTree and computes all the bins
    %   PreAllocateSpace : Pre allocates space for all data for quicker
    %                      computation based on uniform distribution
    %   DeAllocateSpace  : Removes empty pre allocated space
    %   Divide           : Recursive function to divide node
    %   Dividenode       : Computes all variables for given node
    
    properties
        points;
        potentials;
        targets;
        finePoints;
        NN;
        pointIndex;
        nodeCorners;
        nodeLevel;
        nodeParents;
        nodeChildren;
        nodeCount;
        nodeCapacity;
        interactions;
        arguments;
    end
    
    methods
        function this = OcTreeUniform(potPoints,targetPoints,varargin)
            %OcTree Construct an instance of this class
            %   Defines variables for class and initialises the this build
            this.points = unique([potPoints;targetPoints],'rows','stable');
            this.potentials = ismember(this.points,potPoints,'rows');
            this.targets = ismember(this.points,targetPoints,'rows');
            
            numPts = size(this.points,1);
            
            IP = inputParser;
            addOptional(IP,'parThreads',0);
            addOptional(IP,'nodeCapacity',200);
            addOptional(IP,'maxDepth',21);
            addOptional(IP,'finePoints',[]);
            addOptional(IP,'nearest',0);
            addOptional(IP,'NN',[]);
            addOptional(IP,'blockSize',0.2)
            parse(IP,varargin{:});
            this.arguments = IP.Results;
            
            this.nodeCorners = [min(this.points,[],'all')*ones(1,3)...
                               max(this.points,[],'all')*ones(1,3)];
            this.pointIndex = ones(numPts,1);
            this.nodeLevel = 0;
            this.nodeParents(1) = 0;
            this.nodeChildren = 2:9;
            this.nodeCount = 1;
            this.nodeCapacity = min(this.arguments.nodeCapacity,numPts-1);
            
            this.PreAllocateSpace;
            this.Divide;
            this.DeAllocateSpace;
            
            this.geninteractionlists;
            
            if isempty(this.arguments.finePoints)
                this.arguments.nearest = 0;
            else
                this.arguments.nearest = 1;
                this.finePoints = this.arguments.finePoints;
                this.arguments.finePoints = [];
            end
            
            this.NN = [];
            if isempty(this.arguments.NN) && this.arguments.nearest
                this.NN = nearestMatrix(potPoints,this.finePoints,...
                                    this.arguments.blockSize,...
                                    this.arguments.parThreads);
            else
                this.NN = this.arguments.NN;
                this.arguments.NN = [];
            end
        end
        
        function PreAllocateSpace(this)
            % PreAllocateSpace preallocate space in memory for this
            %   Defines empty arrays based on the maximum number of nodes
            %   using the uniform distrubution as the upper limit
            maxLevel = ceil(log(size(this.points,1)...
                       /this.nodeCapacity)/log(8))+1;
            numnodes = sum(8.^(0:maxLevel));
            this.nodeLevel(numnodes) = 0;
            this.nodeParents(numnodes) = 0;
            this.nodeChildren(numnodes,8) = 0;
            this.nodeCorners(numnodes,6) = 0;
        end
        
        function DeAllocateSpace(this)
            % DeAllocateSpace removes extra unused rows from arrays
            %   Sets extra rows in arrays to blank reducing the memory 
            %   used by the this
            
            this.nodeLevel(this.nodeCount+1:end) = [];
            this.nodeParents(this.nodeCount+1:end) = [];
            this.nodeChildren(this.nodeCount+1:end,:) = [];
            this.nodeCorners(this.nodeCount+1:end,:) = [];
        end
        
        function Divide(this)
            %Divide Recursive function to divide points into nodes.
            %   Recursive function which takes an array of nodes and an 
            %   imput and divides them in order. Function will generate 
            %   the this for all nodes below the starting node.
            %
            %Inputs:
            %   startingnodes : Index of nodes which need to be divides.

            this.Dividenode(1);
            level = 1;
            divide = 1;
            nodes = find(this.nodeLevel == level);
            while divide == 1
                divide = 0;
                for i = 1:numel(nodes)
                    this.Dividenode(nodes(i));
                end
                nodes = find(this.nodeLevel == level+1);
                for i = 1:numel(nodes)
                    if nnz(this.pointIndex==nodes(i)) > this.nodeCapacity
                        divide = 1;
                        break;
                    end
                end
                if level+1 >= this.arguments.maxDepth
                    error('Max level reached: Increase Max Depth')
                else
                    level = level + 1;
                end
               
            end
        end
        
        function Dividenode(this,node)
            %Dividenode Recursive fDivide given node and update variables
            %   Divide the given node, computing the children nodes 
            %   indicies, coordinates and update which bin the points are 
            %   in.

            %Inputs:
            %   node : Index of node which need to be divided.
            
            
            nodePtMask = this.pointIndex==node;
            pointsInnode = this.points(nodePtMask,:);
            
            bottomFrontLeft = this.nodeCorners(node,1:3);
            topBackRight = this.nodeCorners(node,4:6);
            
            center = mean([bottomFrontLeft; topBackRight],1);
            coords = [bottomFrontLeft center topBackRight];
            
            newCorners = coords([...
                1 2 3 4 5 6;
                1 2 6 4 5 9;
                1 5 3 4 8 6;
                1 5 6 4 8 9;
                4 2 3 7 5 6;
                4 2 6 7 5 9;
                4 5 3 7 8 6;
                4 5 6 7 8 9]);
            
            nodeMap = cat(3,[0 0 0], [0 0 1], [0 1 0], [0 1 1], [1 0 0], [1 0 1], [1 1 0], [1 1 1]);
            gtMask = pointsInnode > center;
            [~, pointLoc] = max(all(gtMask==nodeMap,2),[],3);
             
            
            newnodeIndex = this.nodeCount+1:this.nodeCount+8;
            this.nodeChildren(node,:) = newnodeIndex;
            this.nodeCorners(newnodeIndex,:) = newCorners;
            this.nodeLevel(newnodeIndex) = this.nodeLevel(node)+1;
            this.nodeParents(newnodeIndex) = node;
            this.pointIndex(nodePtMask) = newnodeIndex(pointLoc);
            this.nodeCount = this.nodeCount+8;
        end
        
        function geninteractionlists(this)
            %GENERATEINTERACTIONLISTS Generate U,V,W and X interaction 
            %                         lists


            this.interactions = cell(this.nodeCount, 4);


            for i=2:this.nodeCount
                [U,V,W]= genUVWinterationlists(this,i);
                this.interactions(i,:)= {U',V',W',[]};
            end


            [M,~] = max(cellfun(@length, this.interactions(:,3)));
            Padded = cellfun(@(x) [x zeros(1, M - numel(x))], ...
                             this.interactions(:,3), 'un', 0);
            W = cat(1,Padded{:});


            for i=2:this.nodeCount
               [X, ~] = find(W==i);
               this.interactions{i,4}= X';
            end

        end
        
        function [out]=getPotentialPoints(this)
            out = this.points(this.potentials,:);
        end
        
        function [out]=getTargetPoints(this)
            out = this.points(this.targets,:);
        end
        
        function [out]=getFinePoints(this)
            out = this.finePoints;
        end
        
    end
    
    methods (Access = private)
       function [U,V,W] = genUVWinterationlists(this, index)
            %GENINTERATIONLIST Generate U,V and W interaction lists for given index
            %Inputs:
            %   index : nodeto compute interaction list for
            %
            %Output:
            %   U : returns a column array off all nodes in U interaction
            %       list
            %   V : returns a column array off all nodes in V interaction
            %       list
            %   W : returns a column array off all nodes in W interaction
            %       list

            indexLevel = this.nodeLevel(index);

            %generate U
            if isleaf(this,index,'Target')
                U = find(nodeintersection(this,index,(1:this.nodeCount)') & isleaf(this,(1:this.nodeCount)','Potential'));
            else
                U = [];
            end

            %generate V
            parentIndex = this.nodeParents(index);
            parentNeigbours = find(nodeintersection(this,parentIndex,(1:this.nodeCount)') & (indexLevel-1 >= this.nodeLevel)');

            parentNeigbours = setdiff(parentNeigbours, [this.nodeParents(parentNeigbours) parentIndex]);

            parentNeigbours = parentNeigbours(~isleaf(this,parentNeigbours,'All'));
            parentNeigboursChildren = reshape(this.nodeChildren(parentNeigbours,:),[],1);
            parentNeigboursChildren = parentNeigboursChildren(~nodeintersection(this,index,parentNeigboursChildren));

            V = setdiff(parentNeigboursChildren,U);


            % generate W
            if isleaf(this,index,'Target')
                neighbours = find(nodeintersection(this,index,(1:this.nodeCount)') & (indexLevel >= this.nodeLevel)');
                neighbours = setdiff(neighbours, this.nodeParents(neighbours));

                neighbourChildren = getChildren(this,neighbours);
                neighbourChildren = neighbourChildren(isleaf(this,neighbourChildren,'Potential'));

                W = setdiff(neighbourChildren, U);
                
                for level = max(this.nodeLevel):-1:1
                    nodes = find(this.nodeLevel==level);
                    for i = 1:size(nodes,2)
                        children = this.nodeChildren(nodes(i),:);
                        if all(ismember(children,W))
                            W = [setdiff(W,children) nodes(i)];
                        end
                    end
                end
                
                
            else
                W = [];
            end
       end
    end
end

