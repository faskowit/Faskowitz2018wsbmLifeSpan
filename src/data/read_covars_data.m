function [covarsDataStruct] = read_covars_data(rawDataDir,datasetDemo,wanted_files)
% read covariate data, provide the project directory, the already read
% dataset so that we can match to the rownames in it, and the names of the
% csv covariate files to read in.
% returns a struct for each of the covar files you want 

docDataDir = strcat(rawDataDir,'/doc/assessment_data/');

disp('reading from the dir:')
disp(docDataDir)

%wanted_files = { '8100_WASI-II_20170606.csv' '8100_WIAT-IIA_20170607.csv' };

covarsDataStruct = struct();

for idx = 1:length(wanted_files)

    readData = readtable(strcat(docDataDir,wanted_files{idx}));
    % remove the first line, which we dont need
    readData(1,:) = [] ;
    % get only subjID unique entries
    [~,uniqIdx] = unique(table2cell(readData(:,1)));
    readData = readData(uniqIdx,:);
    readData.Properties.RowNames = table2cell(readData(:,1)) ;
    
    %readData = standardizeMissing(readData,'');
    %find missing data, remove rows if data is missing (kind of extreme...)
    readData = readData(sum(ismissing(readData),2) == 0,:) ;
    
    [intersectSubj,~,rowIdx] = intersect(readData.Properties.RowNames,...
        datasetDemo.Properties.RowNames);
    
    [~,covarsDataStruct(idx).name] = fileparts(wanted_files{idx});    
    covarsDataStruct(idx).interSubj = intersectSubj ;
    covarsDataStruct(idx).interSubjIdx = rowIdx ;
    covarsDataStruct(idx).interDemo = datasetDemo(intersectSubj,:);
    covarsDataStruct(idx).interCovar = readData(intersectSubj,:);
    
end




