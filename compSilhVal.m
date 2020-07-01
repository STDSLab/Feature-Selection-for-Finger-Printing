function [silhVal, ratio] = compSilhVal(sbjData, dist)
%compSilhVal Compute silhouette value and overlap ratio for a set of vectorized FC scans
%   Description:
%      This function computes the silhouette value and overlap ratio for a set of 
%      vectorized FC scans.
%
%   Input:
%      sbjData - a MATLAB set containing day {1} and {2} vectorized FCs for number of
%       subjects
%      dist - distance string for computing silhouette value (any string that evalclusters
%       supports e.g. 'euclidean'
%
%   Output:
%      silhVal - the average silhouette value for input data
%      ratio - overlap ratio i.e. the ratio of scans with negative silhouette value
%   
%   Author:
%      Kendrick Li [3-29-2019]
    
    nmSbj = size(sbjData{1}, 2);
    E = evalclusters_m([sbjData{1} sbjData{2}]', ...
        [1:nmSbj 1:nmSbj]', 'silhouette', 'Distance', dist);
    subjWthNegScans = unique(E.ObsSilhouetteValues(E.ObsSilhouetteValues(:, 2) < 0));
    ratio = length(subjWthNegScans)/nmSbj;
    silhVal = E.CriterionValues;
end

