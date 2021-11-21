classdef FMM
    %FMM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tree;
        arguments;
    end
    
    properties (Access = private)
        uppot;
        downpot;
    end
    
    methods
        function this = FMM(tree,varargin)
            %FMM Construct an instance of this class
            %   Detailed explanation goes here
            this.tree = tree;
            IP = inputParser;
            addOptional(IP,'epsilon',0.01);
            addOptional(IP,'mu',1);
            addOptional(IP,'GPU',0);
            addOptional(IP,'ParThreads',0);
            addOptional(IP,'BlockSize',1);
            addOptional(IP,'Coronares',6);
            addOptional(IP,'CoronaShells',1);
            parse(IP,varargin{:});
            this.arguments = IP.Results;
        end
        
        function this = vectorProduct(this,potentials)
            
        end
        
    end
    
    methods (Access = private)
        function UpwardPass(this)
            
        end
        
        function DownwardPass(this)
            
        end
        
    end
end

