function [output, sortedScores] = mtfAnalysis(DTn, scores, F, l, transFlatID, flatID, methodNm, dataArgs, tmpData, sessionArgs, analysisArgs)
%mtfAnalysis Perform FC fingerprinting using S edge scores
%   Description:
%      Uses S edge scores to select best nodes and perform FC fingerprinting
%      using the edges of selected nodes only.
%
%   Input:
%      DTn - Number of sets to compare. Should always be equal to
%       2.
%      scores - S edge scores
%      F - the weighted cluster membership matrix
%      l - vector containing different numbers of top S edges to use.
%       set this to a scalar to only perform on one number of top S
%       edges to use. E.g., [3, 6, 9] would return 3 different results
%       using the top 3, top 6, and top 9 S edges.
%      transFlatID - 1 x 2 cell containing the row ID {1} and column ID {2}
%       for each element in the flattened S matrix. Same guidelines as flatID.
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
%      sortedScores - top scores for each S edge where edge scores 
%       not selected are set to 0. Use squareform to see a 2D matrix
%       of selected S edge scores.
%       
%   Author:
%      Kendrick Li [12-12-2019]

    fullData = [tmpData.subjFCs{1} tmpData.subjFCs{2}];
    nmReg = dataArgs.numReg; nmSbj = dataArgs.numSubj;
    nmScans = nmSbj*2; nmEdges = (nmReg*(nmReg - 1))/2;

    fullFC = zeros(nmReg, nmReg, nmScans);
    for iScan = 1:nmScans
      for iEdge = 1:nmEdges
        fullFC(flatID{1}(iEdge), flatID{2}(iEdge), iScan) = ...
          fullData(iEdge, iScan);
      end

      fullFC(:, :, iScan) = ...
        fullFC(:, :, iScan) + fullFC(:, :, iScan)' + eye(nmReg);
    end

    dataArgs.numReg = size(F, 2);
    nmReg = dataArgs.numReg;
    output = zeros(DTn, 3, length(l));
    
    compS = {computeS(fullFC(:, :, 1:nmSbj), F) ...
      computeS(fullFC(:, :, nmSbj+1:end), F)};
    
    flatmap = triu(true(nmReg));
    smlEdges = sum(sum(flatmap));

    flatS = {zeros(smlEdges, nmSbj) ...
      zeros(smlEdges, nmSbj)};

    for iDT = 1:2
      for iScan = 1:nmSbj
        tmpScan = compS{iDT}(:, :, iScan);
        flatS{iDT}(:, iScan) = tmpScan(flatmap);
      end
    end
    
    modData.subjFCs = flatS;
    nodeArgs.flag = true; nodeArgs.methodNm = methodNm; nodeArgs.mode = 2;

    %% do method on test
    % take the top l edges (loop on different sized
    % l) only when doing fingerprinting on test
    % cohort
    [sortedBase, sortedI] = sort(scores);
    for iL = 1:length(l)
        %% create node mask
        nodeArgs.mask = false(nmReg); sorted = sortedBase;
        sorted(l(iL) + 1:end) = 0; sorted(sortedI) = sorted;
        sortedScores = reshape(sorted, size(scores));
        flatMask = sortedScores ~= 0; sel = [transFlatID{1}(flatMask) transFlatID{2}(flatMask)];
        
        for iSel = 1:size(sel, 1)
          nodeArgs.mask(sel(iSel, 1), sel(iSel, 2)) = true;
        end
        
        %% run fingerprinting analysis
        for DT = 1:DTn
            output(DT, :, iL) = maskFP(DT, modData, nodeArgs, ...
                dataArgs, sessionArgs, analysisArgs);
        end
    end
end

