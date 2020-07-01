function [output] = krNodeAnalysis(DTn, ordering, k, r, methodNm, dataArgs, tmpData, sessionArgs, analysisArgs)
%krNodeAnalysis Perform FC fingerprinting using node scores
%   Description:
%      Uses node scores to select best nodes and perform FC fingerprinting
%      using the edges of selected nodes only.
%
%   Input:
%      DTn - Number of sets to compare. Should always be equal to
%       2.
%      ordering - the rank for each node, e.g., if row 5 equals 10,
%       node 5 is the 10th highest scoring node.
%      k - vector containing different numbers of top nodes to use.
%       set this to a scalar to only perform on one number of top
%       nodes to use. E.g., [3, 6, 9] would return 3 different results
%       using the top 3, top 6, and top 9 nodes.
%      r - set to number of regions
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
%       
%   Author:
%      Kendrick Li [11-25-2019]

    nmReg = dataArgs.numReg;
    output = zeros(DTn, 3, length(k));
    
    nodeArgs.flag = true; nodeArgs.methodNm = methodNm; nodeArgs.mode = 3;
    %% use node results for k and r%
    for iK = 1:length(k)
        for DT = 1:DTn
            %% create node mask
            nodeMask = false(nmReg);
            nodeMask(ordering(1:k(iK), DT), :) = true;
            nodeMask(:, ordering(1:k(iK), DT)) = true;

            nodeMask(ordering(r + 1:end, DT), :) = false;
            nodeMask(:, ordering(r + 1:end, DT)) = false;

            nodeMask = triu(nodeMask, 1);
            nodeArgs.mask = nodeMask;

            %% run fingerprinting analysis
            output(DT, :, iK) = maskFP(DT, tmpData, nodeArgs, ...
                dataArgs, sessionArgs, analysisArgs);
        end
    end
end

