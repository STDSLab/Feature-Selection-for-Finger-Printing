classdef   SilhouetteEvaluation_m  < clustering.evaluation.ClusterCriterion
%SilhouetteEvaluation Evaluation based on silhouette criterion
%   MODIFIED to return individual silhouette values.
%   Modification author:
%      Kendrick Li [1-7-2019]
%
%   A SilhouetteEvaluation object computes the silhouette criterion values
%   which can be used for estimating the number of clusters. It's created
%   using the following syntax:
%
%       EVA = EVALCLUSTERS(X,CLUST,'Silhouette')
%
%   SilhouetteEvaluation properties:
%      X                     - Data used to create this object.
%      Distance              - The distance metric used.
%      NumObservations       - Number of observations in X.
%      OptimalK              - The optimal number of clusters suggested.
%      OptimalY              - The optimal clustering solution.
%      InspectedK            - List of the number of clusters inspected.
%      CriterionName         - Name of the criterion.
%      CriterionValues       - The criterion values for each number of clusters.
%      ClusteringFunction    - Clustering algorithm used.
%      Missing               - A vector indicating which observations are missing.
%      ClusterPriors         - The cluster priors.
%      ClusterSilhouettes    - The silhouette values for each cluster in each clustering solution.
%
%   SilhouetteEvaluation methods:
%      addK                  - Add a new list of number of clusters for evaluation.
%      compact               - Compact this object.
%      plot                  - Plot the criterion values.
%
%   See also EVALCLUSTERS, clustering.evaluation.GapEvaluation,
%            clustering.evaluation.CalinskiHarabaszEvaluation,
%            clustering.evaluation.DaviesBouldinEvaluation.

%   References:
%     Kaufman L. and Rousseeuw, P.J. Finding Groups in Data: An
%     Introduction to Cluster Analysis, Wiley, 1990


%  Copyright 2013 The MathWorks, Inc.

    properties(GetAccess=public, SetAccess=protected)
    
%Distance The distance metric used.
%   The Distance property specifies the distance metric used. It's a string
%   representing a built-in distance measure, a function handle, or a
%   numeric distance matrix in the vector form created by PDIST.
        Distance = 'sqEuclidean';
        
%ClusterPriors The cluster priors.
%   The ClusterPriors is one of the strings 'empirical' or 'equal' 
%   indicating the prior probabilities for each class.
        ClusterPriors = 'empirical';
        
%ClusterSilhouettes Silhouette values in each cluster.
%   The ClusterSilhouettes property is a cell array of matrices. The number
%   of cells is equal to the number of clustering solutions. Each cell
%   contains a matrix containing the average silhouette value for each
%   cluster in the corresponding clustering solution.
        ClusterSilhouettes = [];
        
%ObsSilhouetteValues Silhouette values of each observation with a class
% label. Li ADDITION
%   TODO
        ObsSilhouetteValues = [];
  
    end
    
    properties(GetAccess=protected, SetAccess=protected)
        priorLoc;
        %1 means ClusterPriors is 'equal', 2 means ClusterPriors is 'empirical'
    end
    
    methods %Public
               
        function this = addK(this,klist)
            this = addK@clustering.evaluation.ClusterCriterion(this,klist);
          
            [this,newCriterionValues,newClustSilVals] = evalklist(this,...
                this.NewKList);
            
            this.CriterionValues = [this.CriterionValues, newCriterionValues];
            this.ClusterSilhouettes = [this.ClusterSilhouettes, newClustSilVals];
            this.NewKList = [];
        end
    end
    
    methods (Access = protected)
        function [this, criterionValues, clustSilhouettes] = ...
                evalklist(this,klist)
            p = numel(klist);
            criterionValues = nan(1,p);
            clustSilhouettes = cell(1,p);
            if ~isnan(this.OptimalK)
                optimalCritV = max(this.CriterionValues);
            else
                optimalCritV = -inf;
            end
            
            for j = 1:p
                NC = klist(j);
                IDX = this.evalFun(NC);
                if isempty(IDX)
                    continue;
                end
                silVals = silhouette(this.PrivX,IDX,this.Distance);
                
                if ~isa(this.ClusteringFunction,'function_handle')
                    clustSilhouettes{j} = accumarray(IDX, silVals,[NC,1],@mean);
                else
                    clustSilhouettes{j} = accumarray(grp2idx(IDX), silVals,[NC,1],@mean);
                end
                
                if this.priorLoc == 1 %equal
                    criterionValues(j) = mean(clustSilhouettes{j});
                else % empirical prior
                    criterionValues(j) = mean(silVals);
                end
                
                %need to find optimalY,
                if  criterionValues(j) > optimalCritV
                    optimalCritV = criterionValues(j);
                    this.OptimalK = klist(j);
                    this.OptimalY = IDX;
                end
            end
            
            this = postProcess(this);
            
        end %function evalklist
         
    end %protected
    
    methods (Hidden)
        function   this = SilhouetteEvaluation_m(X,Y,varargin)
            narginchk(2,inf);
            [basePnames, baseDflts ] = clustering.evaluation.ClusterCriterion.getAllArgs();
            
            [klist,~,args] = ...
                internal.stats.parseArgs(basePnames, baseDflts, varargin{:});
            
            this = this@clustering.evaluation.ClusterCriterion(X,Y,klist);
            this.CriterionName = 'Silhouette';
             
            % parse input arguments for Silhouette
            pnames = {'Distance'     'ClusterPriors' };
            dflts =  { 'sqeuclidean'  'empirical' };
            
            [distMetric, prior] =...
                internal.stats.parseArgs(pnames, dflts, args{:});
            
            p = numel(this.InspectedK);
            [this.Distance, distLoc] = checkDistance(distMetric, size(this.X,1));
                         
            if isnumeric(this.Distance) % Distance vector
                %Remove the entries with missing values in the distance vector
                if any(this.Missing)
                    this.Distance = squareform(this.Distance,'tomatrix');
                    this.Distance=this.Distance(~this.Missing,~this.Missing);
                    this.Distance = squareform(this.Distance,'tovector');      
                end
                if any(isnan(this.Distance))
                    error(message('stats:clustering:evaluation:SilhouetteEvaluation:TooManyNaNsInDist'));
                end
            end
            
            if iscell(Y)%if Y is soft clustering solutions, convert it hard solutions
               Y = soft2hard(this,Y);
            end
          
            %Sanity check for clusterPrior
            [this.ClusterPriors,this.priorLoc]=internal.stats.getParamVal(...
                prior,{'equal','empirical'},'ClusterPrior');
                                  
            %loop over on each clustering solution
            if isnumeric(Y) % Y is a matrix
                Y(this.Missing,:) = [];
                this.ClusterSilhouettes = cell(1,p);
                for j = 1:p
                    NC = this.InspectedK(j);
                    %Apply the specified distance metric.
                    silVals = silhouette(this.PrivX, Y(:,j), this.Distance);
                    this.ClusterSilhouettes{j} = accumarray(grp2idx(Y(:,j)), silVals,[NC,1],@mean);
                    this.ObsSilhouetteValues = [Y silVals]; % Li ADDITION
                    if this.priorLoc == 1 %equal
                        % Each cluster makes the same contribution
                        this.CriterionValues(j) = mean(this.ClusterSilhouettes{j});
                    else %empirical
                        this.CriterionValues(j) = mean(silVals);
                    end 
                end
              
                [maxVal,optimalKLoc] = max(this.CriterionValues);
                if ~isnan(maxVal)
                    this.OptimalK = this.InspectedK(optimalKLoc);
                end
                
                this.PrivX = [];
                
            else % Y represents a clustering algorithm
                if ischar(Y) %Y is a built-in function
                    if distLoc ~= 1 % distance is not squared Euclidean distance
                        %Apply the specified distance metrics on built-in functions
                        if this.FunLoc == 1 %kmeans
                             this.ClustFuncHandle = ...
                                @(X,NC)(kmeans(X,NC,'rep',5,'distance',this.Distance,'emptyaction','singleton'));
                        elseif this.FunLoc == 2 %linkage
                            %For linkage, distance must be one of the string in
                            %distList or a function handle or a distance matrix
                            if isnumeric(this.Distance) %distance is a distance matrix
                                %Use average linkage since we don't know
                                %whether the distance comes from the Euclidean
                                %distance and we try to avoid calling expensive
                                %"isEuclidean"
                               this.LinkageFuncHandle = [];
                               this.LinkageZ = linkage(this.Distance,'ave');
                                
                            elseif distLoc ~= 2% not 'euclidean'
                                %the distance is neither 'euclidean' nor 'sqeuclidean'
                                this.LinkageFuncHandle = @(X)linkage(X,'ave',this.Distance);
                            end
                        elseif this.FunLoc == 3 %gmdistribution
                            %For gmdistribution, it must be sqEuclidean distance
                            error(message('stats:clustering:evaluation:SilhouetteEvaluation:BadDistForGM'));
                        end
                    end % distLoc ~= 1
                    
                    if  this.FunLoc == 2 && isempty(this.LinkageZ) %linkage with squared EucLidean distance
                        this.LinkageZ = this.LinkageFuncHandle(this.PrivX);
                    end
                end
                
                [this,this.CriterionValues,this.ClusterSilhouettes] = ...
                    evalklist(this, this.InspectedK);
            end
        end %constructor
       
  end
end

function [distance,distLoc] = checkDistance(distance,nIn)
if ischar(distance)
    distList = {'sqEuclidean' 'Euclidean' 'cityblock'  'cosine' 'correlation' 'Hamming', 'Jaccard'};
   [distance,distLoc] = internal.stats.getParamVal(distance,distList,'distance');
elseif isnumeric(distance)
    
    if (size(distance,1) == 1) && (size(distance,2) == .5*nIn*(nIn-1))
        if any(distance < 0)
            error(message('stats:clustering:evaluation:SilhouetteEvaluation:NegDistanceValues'));
        end
    else
        error(message('stats:silhouette:DistanceMatrixNotUpperTri'));
    end
    distLoc = 0;
elseif isa(distance,'function_handle')
       distLoc = 0;
else
    error(message('stats:clustering:evaluation:SilhouetteEvaluation:BadDist'));
end
end
