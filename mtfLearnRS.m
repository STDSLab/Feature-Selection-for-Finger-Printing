function [F, scores] = mtfLearnRS(flatS, nmTrainSbj, F, smlEdges, nmScans)
%mtfLearnRS Scores edges in S based on rank sum cost
%   Description:
%      Scores each edge in S based on the rank sum (RS) cost.
%
%   Input:
%      flatS - a 1 x 2 cell where each contains a smlEdges x nmSbj matrix 
%       containing vectorized S matrices.
%      nmSbj - number of subjects
%      F - the weighted cluster membership matrix - unused
%      smlEdges - number of edges used from the S matrix, i.e., number of
%       unique values in the S matrix
%      nmScans - number of scans, should be equal to nmSbj*2
%
%   Output:
%      F - returns F unchanged
%      scores - the score for each unique edge in S
%       
%   Author:
%      Kendrick Li [12-10-2019]

    interlFlatFC = flatS{2}(:, [1; 1]*(1:nmTrainSbj));
    interlFlatFC(:, 1:2:end) = flatS{1};
    interlFlatFC = interlFlatFC.';
    
    %% find best edges by calculating the rank sum
    scores = zeros(smlEdges, 1);
    for iEdge = 1:smlEdges
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
    end
end