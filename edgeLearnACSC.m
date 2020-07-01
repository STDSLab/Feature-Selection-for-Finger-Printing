function [scores] = edgeLearnACSC(nmEdges, nmSbj, nmScans, flatFCdata)
%edgeLearnACSC Scores edges based on average cross-session cost
%   Description:
%      Scores each edge based on the average cross-session cost (ACSC).
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

    diagSelect = diag(true(nmSbj, 1));

    %% find best edges by calculating the average cross-session cost
    scores = zeros(nmEdges, 1);
    for iEdge = 1:nmEdges
      simMat = pdist2(flatFCdata{1}(iEdge, :)', ...
        flatFCdata{2}(iEdge, :)', 'seuclidean');

      % calc new cost
      scores(iEdge) = mean(simMat(diagSelect)) - ...
        mean(simMat(~diagSelect));
    end
end

