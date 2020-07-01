function [scores] = edgeLearnRS(nmEdges, nmSbj, nmScans, flatFCdata)
%edgeLearnRS Scores edges based on rank sum cost
%   Description:
%      Scores each edge based on the rank sum (RS) cost.
%
%   Input:
%      nmEdges - number of edges. Computed by (N*(N-1))/2
%      nmSbj - number of subjects
%      nmScans - number of scans, should be equal to nmSbj*2
%      flatFCdata - a 1 x 2 cell where each contains a nmEdges x nmSbj 
%       matrix containing unique edge values. Each column is one scan.
%
%   Output:
%      scores - edge scores
%       
%   Author:
%      Kendrick Li [12-2-2019]

    interlFlatFC = flatFCdata{2}(:, [1; 1]*(1:nmSbj));
    interlFlatFC(:, 1:2:end) = flatFCdata{1};
    interlFlatFC = interlFlatFC.';
    
    %appFlatFC = [flatFCdata{1} flatFCdata{2}].';
    
    clear flatFCdata;

    %% find best edges by calculating the rank sum
    scores = zeros(nmEdges, 1);
    for iEdge = 1:nmEdges
        %% calc edge dist
        distMat = squareform(pdist(interlFlatFC(:, iEdge), 'seuclidean'));

        %% generate ranks
        for iRMCol = 2:2:nmScans
            intraSubjDist = distMat(iRMCol - 1, iRMCol);

            %% rank sum
            scores(iEdge) = scores(iEdge) + ...
                sum(distMat(iRMCol - 1, :) < intraSubjDist);
            scores(iEdge) = scores(iEdge) + ...
                sum(distMat(iRMCol, :) < intraSubjDist);
        end
        
        %{
        distMat2 = squareform(pdist(appFlatFC(:, iEdge), 'seuclidean'));
        rankMat = zeros(size(distMat2));
        for iR = 1:nmScans
          [~, order] = sort(distMat2(iR, :));
          rankMat(iR, order) = 0:nmScans-1;
        end
        
        secondERS = sum(diag(rankMat(11:20, 1:10)));
        secondERS = secondERS + sum(diag(rankMat(1:10, 11:20)));
        if scores ~= secondERS
          warning('VALUES DO NOT MATCH');
        end
        %}
    end
end