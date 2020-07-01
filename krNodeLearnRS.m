function [scores, ordering] = krNodeLearnRS(DTn, dataArgs, flatID, flatFCdata, distanceMetric, sessionArgs, analysisArgs)
%krNodeLearnRS Scores and orders nodes based on rank sum cost
%   Description:
%      Scores each node based on the rank sum (RS) cost
%      and orders nodes based on the best scores.
%
%   Input:
%      DTn - Number of sets to compare. Should always be equal to
%       2.
%      dataArgs - data arguments, contains:
%       numSubj - number of subjects
%       numReg - number of regions (nodes)
%      flatID - 1 x 2 cell containing the row ID {1} and column ID {2}
%       for each element in the flattened FC. E.g. for a 3x3 FC matrix,
%       if the FC was flattened column-wise, flatID{1} = [1, 1, 2] and 
%       flatID{2} = [2, 3, 3] (corresponding to edges [1, 2], [1, 3], 
%       and [2, 3]). If the diagonal is included, flatID should also 
%       include IDs for the diagonal edges.
%      flatFCdata - a 1 x 2 cell where each contains a nmEdges x nmSbj 
%       matrix containing unique edge values. Each column is one scan.
%       Make sure that flatID matches this matrix, i.e., if FCs were 
%       flattened column-wise flatID should be ordered column-wise.
%      distanceMetric - distance metric string for pdist2. Any distance
%       string that pdist2 supports can be used, e.g., 'correlation'.
%      sessionArgs - unused
%      analysisArgs - unused
%
%   Output:
%      scores - ACSC scores for each node. As ACSC is identical regardless
%       of whether the first set is used as the database or the second set,
%       both columns are identical.
%      ordering - contains the rank for each node, e.g., if row 5 equals 10,
%       node 5 is the 10th highest scoring node.
%       
%   Author:
%      Kendrick Li [12-3-2019]

    nmReg = dataArgs.numReg; nmSbj = dataArgs.numSubj; nmScans = nmSbj*2;

    interlFlatFC = flatFCdata{2}(:, [1; 1]*(1:nmSbj));
    interlFlatFC(:, 1:2:end) = flatFCdata{1};
    interlFlatFC = interlFlatFC.';
    
    clear flatFCdata;
    
    %% do training single node fingerprinting
    scores = zeros(nmReg, DTn);
    for iNode = 1:nmReg
      % extract node data from flatFC
      combMask = flatID{1} == iNode | flatID{2} == iNode;

      % gen new subject session distance
      distMat = squareform(pdist(interlFlatFC(:, combMask), distanceMetric));

      % calc new cost
      for iRMCol = 2:2:nmScans
        intraSubjDist = distMat(iRMCol - 1, iRMCol);

        %% rank sum
        scores(iNode, 1) = scores(iNode, 1) + ...
          sum(distMat(iRMCol - 1, :) < intraSubjDist);
        scores(iNode, 1) = scores(iNode, 1) + ...
          sum(distMat(iRMCol, :) < intraSubjDist);
      end
      scores(iNode, 1) = -scores(iNode, 1);
      scores(iNode, 2) = scores(iNode, 1);
    end

    %% order node results
    ordering = zeros(nmReg, DTn);
    [~, ordering(:, 1)] = sort(scores(:, 1), 'descend');
    ordering(:, 2) = ordering(:, 1);
    %[~, ordering(:, 2)] = sort(scores(:, 2), 'descend');
end

