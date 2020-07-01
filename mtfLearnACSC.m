function [F, scores] = mtfLearnACSC(flatS, nmSbj, F, smlEdges)
%mtfLearnACSC Scores edges in S based on average cross-session cost
%   Description:
%      Scores each edge in S based on the average cross-session cost (ACSC).
%
%   Input:
%      flatS - a 1 x 2 cell where each contains a smlEdges x nmSbj matrix 
%       containing vectorized S matrices.
%      nmSbj - number of subjects
%      F - the weighted cluster membership matrix - unused
%      smlEdges - number of edges used from the S matrix, i.e., number of
%       unique values in the S matrix
%
%   Output:
%      F - returns F unchanged
%      scores - the score for each unique edge in S
%       
%   Author:
%      Kendrick Li [12-2-2019]
    
    diagSelect = diag(true(nmSbj, 1));

    %% find best edges by calculating the average cross-session cost
    scores = zeros(smlEdges, 1);
    for iEdge = 1:smlEdges
      simMat = pdist2(flatS{1}(iEdge, :)', ...
        flatS{2}(iEdge, :)', 'seuclidean');

      % calc new cost
      scores(iEdge) = mean(simMat(diagSelect)) - ...
        mean(simMat(~diagSelect));
    end
end

