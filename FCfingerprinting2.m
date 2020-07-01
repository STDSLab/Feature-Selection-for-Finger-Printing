function [ output ] = FCfingerprinting2( dataArgs, sessionArgs, analysisArgs )
%FCfingerprinting Run FC fingerprinting on specified data input
%   Description:
%      This function implements FC fingerprinting based on Finn et al. 2015. It performs
%      returns FC fingerprinting on the input data and returns the accuracy and results.
%
%      Finn, Emily S., et al. "Functional connectome fingerprinting: identifying 
%       individuals using patterns of brain connectivity." Nature neuroscience 18.11 
%       (2015): 1664-1671.
%
%   Input:
%      dataArgs - data arguments, contains:
%	    numSubj - number of subjects
%	    numReg - number of regions, i.e., number of nodes
%	    subjData - a 1 x 2 cell where each contains flattened selected FC scans.
%        Rows are selected edges and columns are number of subjects. {1} contains
%        the day 1 scan set and {2} contains the day 2 scan set.
%      sessionArgs - session arguments, contains:
%       database - which set to use as the database. The other scan set will be used
%        as target scans. Can be set to 1 (use first day) or 2 (use second day)
%      analysisArgs - analysis arguments, contains:
%	    flag - whether to do analysis or not, by default set to 0
%
%   Output:
%      output - an output object which contains:
%       dataTargetComparisonMat - the pairwise correlation matrix
%       identification - identification for all target scans
%       identificationTargetAsDatabase - identification using target
%        scans as database and database scans as target
%       idAcc - identification accuracy using (1) database as database
%        and (2) using database as target
%       database - what was used as the database
%       target - what was used as the target scans
%   
%   Author:
%      Kendrick Li [8-5-2019]

    %% HCP data
    numReg = dataArgs.numReg; numSubj = dataArgs.numSubj;
    
    % if we need to use a subset of the entire data, take a random sampling
    % of the dataset
    %{
    if dataArgs.forceNumSubj && dataArgs.numSubj < numSubj
        randomOrder = randperm(numSubj);
        fullData.subjFileNames = fullData.subjFileNames(randomOrder);
        fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeries(randomOrder, :);
        numSubj = dataArgs.numSubj;
        fullData.subjFileNames = fullData.subjFileNames(1:numSubj);
        fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeries(1:numSubj, :);
        if dataArgs.preLoadedFlag
            fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeriesConcat(randomOrder, :);
            fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeriesConcat(1:numSubj, :);
        end
        if dataArgs.preFCFlag
            fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeriesFCAvg(randomOrder, :);
            fullData.subjRegionTimeSeries = fullData.subjRegionTimeSeriesFCAvg(1:numSubj, :);
        end
    end
    %}
    
    % half of the targets are used as the database, so switch the target
    % scans to the database locations
    % in effect removes sampling day bias
    %{
    if dataArgs.randomizeDataTarg
        randomOrder = randperm(numSubj);
        randomOrder(randomOrder > numSubj/2) = 0;
        output.randomOrder = randomOrder;
        for i = 1:numSubj
            if randomOrder(i) == 0
                tmp = cell(1, 4);
                tmp(1) = fullData.subjRegionTimeSeries(i, 3);
                tmp(2) = fullData.subjRegionTimeSeries(i, 4);
                tmp(3) = fullData.subjRegionTimeSeries(i, 1);
                tmp(4) = fullData.subjRegionTimeSeries(i, 2);
                
                fullData.subjRegionTimeSeries(i, :) = tmp;
            end
        end
    end
    %}
       
    % specify to use which rest day (rest 1 or rest 2)
    if sessionArgs.database == 1
        databaseSession = 1;
        targetSession = 2;
    else
        databaseSession = 2;
        targetSession = 1;
    end
    
    % if square matrix we won't be able to tell which is which
    if size(dataArgs.subjData{databaseSession}, 1) == numSubj && ...
       size(dataArgs.subjData{databaseSession}, 2) ~= numSubj
        FCdatabase = dataArgs.subjData{databaseSession}';
        FCtarget = dataArgs.subjData{targetSession}';
    else
        FCdatabase = dataArgs.subjData{databaseSession};
        FCtarget = dataArgs.subjData{targetSession};
    end
    
    %% Fingerprinting
    % calculate the correlation between every target with every database point
    dataTargetComparisonMat = corr(FCdatabase, FCtarget);

    % take the max correlation of the columns as each column is a target, we
    % compare one target to all database values and then choose the database
    % index where the correlation is highest
    [~, identification] = max(dataTargetComparisonMat);

    % if we reverse that and instead take the max correlation of the rows, we
    % are comparing each database to every target, thus reversing the previous 
    % method and using the targets as the database instead
    [~, identificationTargetAsDatabase] = max(dataTargetComparisonMat, [], 2);

    idAcc = [sum(identification == 1:numSubj)/numSubj ...
        sum(identificationTargetAsDatabase' == 1:numSubj)/numSubj];

    output.dataTargetComparisonMat = dataTargetComparisonMat;
    output.identification = identification;
    output.identificationTargetAsDatabase = identificationTargetAsDatabase';
    output.idAcc = idAcc;
    
    %% analysis
    if analysisArgs.flag
        if analysisArgs.method == 1
            output.std = std(FCdatabase, 0, 1);
            output.diff = FCdatabase' - FCtarget';
        elseif analysisArgs.method == 2
            
        end
    end
    
    output.database = FCdatabase'; output.target = FCtarget';
end

