classdef OcTree < handle
    %OCTREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        points;
        pointIndex;
        binCorners;
        binLevel;
        binParents;
        binChildren;
        binCount;
        arguments;
    end
    
    methods
        function this = OcTree(pts,varargin)
            numPts = size(pts,1);
            
            IP = inputParser;
            addOptional(IP,'binCapacity',ceil(numPts/100));
            addOptional(IP,'maxDepth',21);
            parse(IP,varargin{:});
            this.arguments = IP.Results;
            
            this.binCorners = [min(pts,[],1) max(pts,[],1)];
            this.points = pts;
            this.pointIndex = ones(numPts,1);
            this.binLevel = 0;
            this.binParents(1) = 0;
            this.binChildren = 2:9;
            this.binCount = 1;
            
            this.PreAllocateSpace;
            this.Divide(1);
            this.DeAllocateSpace;
        end
        
        function PreAllocateSpace(this)
            numBins = ceil(8*size(this.points,1)/this.arguments.binCapacity);
            this.binLevel(numBins) = 0;
            this.binParents(numBins) = 0;
            this.binChildren(numBins,8) = 0;
            this.binCorners(numBins,6) = 0;
        end
        
        function DeAllocateSpace(this)
            this.binLevel(this.binCount+1:end) = [];
            this.binParents(this.binCount+1:end) = [];
            this.binChildren(this.binCount+1:end,:) = [];
            this.binCorners(this.binCount+1:end,:) = [];
        end
        
        function Divide(this, startingBins)
            for i = 1:length(startingBins)
               bin = startingBins(i);
               
               if this.binLevel(bin)+1 >= this.arguments.maxDepth
                   error('Max level reached: Increase Max Depth')
               end
               
               oldCount = this.binCount;
               if nnz(this.pointIndex==bin) > this.arguments.binCapacity
                   this.DivideBin(bin);
                   this.Divide(oldCount+1:this.binCount);
                   continue;
               end
            end
        end
        
        function DivideBin(this,bin)
            % get points in bin
            binPtMask = this.pointIndex==bin;
            pointsInBin = this.points(binPtMask,:);
            
            bottomFrontLeft = this.binCorners(bin,1:3);
            topBackRight = this.binCorners(bin,4:6);
            
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
            
            binMap = cat(3,[0 0 0], [0 0 1], [0 1 0], [0 1 1], [1 0 0], [1 0 1], [1 1 0], [1 1 1]);
            gtMask = pointsInBin > center;
            [~, pointLoc] = max(all(gtMask==binMap,2),[],3);
             
            
            newBinIndex = this.binCount+1:this.binCount+8;
            this.binChildren(bin,:) = newBinIndex;
            this.binCorners(newBinIndex,:) = newCorners;
            this.binLevel(newBinIndex) = this.binLevel(bin)+1;
            this.binParents(newBinIndex) = bin;
            this.pointIndex(binPtMask) = newBinIndex(pointLoc);
            this.binCount = this.binCount+8;
        end
    end
end

