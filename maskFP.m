function [output] = maskFP(DT, tmpData, nodeArgs, dataArgs, sessionArgs, analysisArgs)
%maskFP Perform masked FC fingerprinting on input data
%   Description:
%      Perform FC fingerprinting using only a subset of edges in the FC.
%
%   Input:
%      DT - which session to work on
%      tmpData - object which contains:
%       subjFCs - a 1 x 2 cell where each contains a nmEdges x nmSbj 
%        matrix containing unique edge values. Each column is one scan.
%      nodeArgs - arguments for selecting which edges to use in the FC
%      dataArgs - see ~Analysis functions
%      sessionArgs - session arguments, see FCfingerprinting2
%      analysisArgs - analysis arguments, see FCfingerprinting2
%
%   Output:
%      output - output vector containing the accuracy, silhouette value, 
%       and overlap ratio
%       
%   Author:
%      Kendrick Li [12-20-2019]

    output = zeros(1, 3);
    
    %% run fingerprinting
    dataArgs.subjData = ...
        {genFCpnts2(tmpData.subjFCs(:, 1), dataArgs.numReg, nodeArgs) ...
         genFCpnts2(tmpData.subjFCs(:, 2), dataArgs.numReg, nodeArgs)};
     
    if numel(dataArgs.subjData{1}) ~= 0
        FPres = FCfingerprinting2(dataArgs, sessionArgs, analysisArgs);
    else
        FPres.idAcc = [-1, -1];
    end
    output(1) = FPres.idAcc(DT);

    %% do silh val change
    if numel(dataArgs.subjData{1}) ~= 0
        %[output(2), output(3)] = compSilhVal(dataArgs.subjData, ...
        %    'correlation');
        [output(2), output(3)] = compSilhVal(dataArgs.subjData, ...
            'euclidean');
    else
        output(3) = -1; output(2) = -1;
    end
end

