function [ FCpnts ] = flattenMatrixToVectorFC( FC, numRegions, mode )
%flattenMatrixToVectorFC Converts FC matrix into an FC vector
%   Description:
%      Vectorizes the input FC matrix using a mask. Can work on any matrix.
%
%   Input:
%      FC - N x N matrix to flatten
%      numRegions - N, the total number of regions (i.e. number of rows),
%       unused
%      mode - 'optional', specifies flattening mode:
%       1 - column-wise flattening, excluding diagonal
%       2 - row-wise flattening, including diagonal
%       3 - row-wise flattening, excluding diagonal (default)
%       4 - no flattening (output = input)
%
%   Output:
%      FCpnts - a vectorized FC based on the input options
%       
%   Author:
%      Kendrick Li [12-13-2019]

    %FCpnts = zeros(numRegions*(numRegions - 1)/2, 1);                   % regionCovariances
    squareFlag = true;
    if nargin == 3 && mode == 1
        tMask = tril(true(size(FC)), -1);
    elseif nargin == 3 && mode == 2
        tMask = triu(true(size(FC)));
    elseif nargin == 3 && mode == 4
        squareFlag = false;
    else
        tMask = triu(true(size(FC)), 1);
    end
    
    if squareFlag
        FCpnts = FC(tMask);
    else
        FCpnts = FC;
    end
      
    
    %{
    cursor = 1;
    for region = 1:numRegions - 1                                               % for each region in matrix (except last)
        nextCursor = cursor - 1 + numRegions - region;
        FCpnts(cursor:nextCursor) = FC(region, region + 1:numRegions); % append unique covariances to the covariance list
        cursor = nextCursor + 1;
    end
    %}
end

