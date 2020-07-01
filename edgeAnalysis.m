function [output, sortedRanks] = edgeAnalysis(DTn, scores, l, flatID, methodNm, dataArgs, tmpData, sessionArgs, analysisArgs)
%edgeAnalysis Perform FC fingerprinting using edge scores
%   Description:
%      Uses edge scores to select best nodes and perform FC fingerprinting
%      using the edges of selected nodes only.
%
%   Input:
%      DTn - Number of sets to compare. Should always be equal to
%       2.
%      scores - edge scores
%      l - vector containing different numbers of top edges to use.
%       set this to a scalar to only perform on one number of top
%       edges to use. E.g., [3, 6, 9] would return 3 different results
%       using the top 3, top 6, and top 9 edges.
%      flatID - 1 x 2 cell containing the row ID {1} and column ID {2}
%       for each element in the flattened FC. E.g. for a 3x3 FC matrix,
%       if the FC was flattened column-wise, flatID{1} = [1, 1, 2] and 
%       flatID{2} = [2, 3, 3] (corresponding to edges [1, 2], [1, 3], 
%       and [2, 3]). If the diagonal is included, flatID should also 
%       include IDs for the diagonal edges.
%      methodNm - masking mode for FC flattening, see genFCpnts2
%      dataArgs - data arguments, contains:
%       numSubj - number of subjects
%       numReg - number of regions, i.e. number of nodes
%      tmpData - object which contains:
%       subjFCs - a 1 x 2 cell where each contains a nmEdges x nmSbj 
%        matrix containing unique edge values. Each column is one scan.
%      sessionArgs - session arguments, see FCfingerprinting2
%      analysisArgs - analysis arguments, see FCfingerprinting2
%
%   Output:
%      output - a DTn x 3 x length(k) matrix where output(i, :, j)
%       contains the accuracy, silhouette value, and overlap ratio
%       for session i (either day 1 database or day 2 database) and
%       each number of top nodes in k.
%      sortedRanks - top scores for each edge where edge scores not
%       selected are set to 0. Use squareform to see a 2D matrix
%       of selected edge scores.
%       
%   Author:
%      Kendrick Li [12-2-2019]

    nmReg = dataArgs.numReg;
    output = zeros(DTn, 3, length(l));
    
    nodeArgs.flag = true; nodeArgs.methodNm = methodNm; nodeArgs.mode = 3;

    %% do method on test
    % take the top l edges (loop on different sized
    % l) only when doing fingerprinting on test
    % cohort
    [sortedBase, sortedI] = sort(scores);
    for iL = 1:length(l)
        %% create node mask
        nodeArgs.mask = false(nmReg); sorted = sortedBase;
        sorted(l(iL) + 1:end) = 0; sorted(sortedI) = sorted;
        sortedRanks = reshape(sorted, size(scores));
        flatMask = sortedRanks ~= 0; sel = [flatID{1}(flatMask) flatID{2}(flatMask)];
        nodeArgs.mask((sel(:, 2) - 1)*nmReg + sel(:, 1)) = true;
        
        %% run fingerprinting analysis
        for DT = 1:DTn
            output(DT, :, iL) = maskFP(DT, tmpData, nodeArgs, ...
                dataArgs, sessionArgs, analysisArgs);
        end
    end
end

