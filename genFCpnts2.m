function [ FCs ] = genFCpnts2( FCscans, numReg, nodeArgs )
%genFCpnts [ out ]. For an input population FCs, flatten FCs
%   Description:
%      Vectorizes all input FC matrices and selects which FC
%      elements to use based on input arguments.
%
%   Input:
%      FCscans - 1 x 1 cell containing flattened matrices (edges x scan)
%      numReg - N, the total number of regions (i.e. number of rows)
%      nodeArgs - node arguments, contains:
%       flag - whether to do special flatten or selection of edges
%       methodNm - method name string, can be:
%        'unsorted' - depreciated
%        'sorted' - depreciated
%        'mask' - use 2D mask
%        'maskF' - use flat mask
%       mode - mask mode, see flattenMatrixToVectorFC
%       mask - the 2D mask which selects which elements to use in the FC
%
%   Output:
%      FCs - vectorized FC matrix based on the input options
%       
%   Author:
%      Kendrick Li [11-25-2019]

    if nargin <= 2
        nodeFlag = false;
    else
        nodeFlag = nodeArgs.flag;
        if nodeFlag
            methodNm = nodeArgs.methodNm;
            
            if strcmp(methodNm, 'unsorted')
                method = 1;
                node = nodeArgs.node;
            elseif strcmp(methodNm, 'sorted')
                method = 2;
                k = nodeArgs.k;
                r = nodeArgs.r;
                
                if k >= r
                    disp('K is bigger than or equal to r ERROR');
                end
            elseif strcmp(methodNm, 'mask')
                method = 3;
                mask = nodeArgs.mask;
            elseif strcmp(methodNm, 'maskF')
                method = 4;
                %% ADDED CHANGE TO INCLUDE FLATTEN MODE
                %   Will likely break all prior experiments
                %   If you are crashing here, add nodeArgs.mode = 1;
                %    to your nodeArgs to make everything work again
                mask = flattenMatrixToVectorFC(nodeArgs.mask, ...
                    numReg, nodeArgs.mode);
            else
                method = 0;
                disp('Unknown method');
            end
        end
    end

    % numScan will always be the whole population
    numScan = length(FCscans);
    if ~nodeFlag
        FCs = zeros(numScan, numReg*(numReg - 1)/2);

        % for each scan, calculate the region correlations and flatten them
        % into vectors (use squareform to generate FCs again)
        for scan = 1:numScan
            FCs(scan, :) = flattenMatrixToVectorFC(FCscans{scan}, numReg);
        end
    else
        if method == 1
            FCs = zeros(numScan, numReg - 1);

            for scan = 1:numScan
                FC = FCscans{scan};
                
                FCs(scan, :) = FC(node, FC(node, :) < 1);
                
                %{
                nodeVals = FC(node, :);
                nodeVals(node) = [];
                FCs(scan, :) = nodeVals;
                %}
            end
        elseif method == 2
            leftk = numReg - k; leftr = numReg - r;
            edges = numReg*(numReg - 1)/2 - leftk*(leftk - 1)/2;
            edges = edges - leftr*k;
            FCs = zeros(numScan, edges);
            
            krMask = true(numReg);
            krMask(k + 1:end, :) = false;
            krMask(:, r + 1:end) = false;
            krMask = triu(krMask, 1);
            
            for scan = 1:numScan
                FC = FCscans{scan};
                
                %{
                cursor = 1;
                for i = 1:k
                    cursorend = cursor - 1 + numReg - i - leftr;
                    FCs(scan, cursor:cursorend) = FC(i, i + 1:r);
                    cursor = cursorend + 1;
                end
                %}
                
                FCs(scan, :) = FC(krMask);
            end
        elseif method == 3
            FCs = zeros(numScan, sum(sum(mask)));
            
            for scan = 1:numScan
                FCs(scan, :) = FCscans{scan}(mask);
            end
        elseif method == 4
            flatFCscans = FCscans{1};
            FCs = flatFCscans(mask, :);
        else
            FCs = null;
        end
    end
end

