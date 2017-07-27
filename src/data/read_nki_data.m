function datasetStruct = read_nki_data(rawDataDir , parcellation )
% read the data from raw_data directoy and put it into datasetStruct 
%
% this function will assume that the data is organized in a specific manner
% that can be read in by this function

%% 
disp('reading from the dir:')
disp(rawDataDir)

%% demographics 

demographicsRaw = tdfread(strcat(rawDataDir,'/participants.tsv')) ;
demographicsRaw = struct2table(demographicsRaw) ; 
demographicsRaw.Properties.VariableNames{'participant_id'} = 'ids' ;
demographicsRaw.Properties.RowNames = cellstr(demographicsRaw.ids) ;

%% connectome data

dataRaw = struct() ;

if strcmp(parcellation, 'yeo')
    parcelGlob = 'yeo' ; 
elseif strcmp(parcellation, 'fs150')
    parcelGlob = 'fs150';
elseif strcmp(parcellation, 'scale125')
    parcelGlob = 'scale125';
else
    disp('invalide parcellation string')
    return 
end

disp(parcelGlob)

% dirs
countMatsDir = strcat(rawDataDir, '/', parcelGlob, '/counts/')  ;
volMatsDir = strcat(rawDataDir, '/', parcelGlob, '/vols/') ;
coorMatsDir = strcat(rawDataDir, '/', parcelGlob, '/coordsMM/') ;
lensMatsDir = strcat(rawDataDir, '/', parcelGlob, '/lens/') ;
    
% record which subjects have complete data...
% we will only deal with them 
ids = cellstr(demographicsRaw.ids) ; 
filterIncompleteData = zeros(length(ids),1) ;
for idx=1:length(ids) 
    
    disp(idx) ;

    tempPotFile = dir(char(strcat(countMatsDir,'*',ids(idx),'*.csv'))) ;    
    if exist(strcat(countMatsDir,tempPotFile.name),'file' ) == 2

        % now look for the vol file as well
        tempPotFile = dir(char(strcat(volMatsDir,'*',ids(idx),'*.csv'))) ;    
        if exist(strcat(volMatsDir,tempPotFile.name),'file' ) == 2
            
            % now look for the coor file 
            tempPotFile = dir(char(strcat(coorMatsDir,'*',ids(idx),'*.csv'))) ;
            if exist(strcat(coorMatsDir,tempPotFile.name),'file' ) == 2
            
                % now look for the coor file 
                tempPotFile = dir(char(strcat(lensMatsDir,'*',ids(idx),'*.csv'))) ;
                if exist(strcat(lensMatsDir,tempPotFile.name),'file' ) == 2

                    disp('found all')
                    filterIncompleteData(idx) = 1 ;
                end
            end
        end
    end
end

% filter out incomplete data
demographicsRaw = demographicsRaw(filterIncompleteData == 1,:) ; 
ids = cellstr(demographicsRaw.ids) ; 

%% actually read in data

for idx=1:length(ids)
          
    % add id colum 
    dataRaw(idx).id = ids(idx) ;
    
    disp(idx)
    
    %% read in counts
    %%%%%%%%%%%%%%%%%
    
    tempPotFile = dir(char(strcat(countMatsDir,'*',ids(idx),'*.csv'))) ; 
    % put into struct
    fileName = char(strcat(countMatsDir,tempPotFile.name));
    % read in file
    dataRaw(idx).countMat = csvread(fileName) ; 
    
    %counts mat
    % get rid of row+col 1
    % because dipy records the connection to '0'
    dataRaw(idx).countMat(1,:) = [] ;
    dataRaw(idx).countMat(:,1) = [] ;
    
    %% read in lengths
    %%%%%%%%%%%%%%%%%%
    
    tempPotFile = dir(char(strcat(lensMatsDir,'*',ids(idx),'*.csv'))) ; 
    % put into struct
    fileName = char(strcat(lensMatsDir,tempPotFile.name));
    % read in file
    dataRaw(idx).lensMat = csvread(fileName) ; 
    
    % fix lens mat
    temp = triu(dataRaw(idx).lensMat,1) ;
    dataRaw(idx).lensMat = dataRaw(idx).lensMat + temp' ;
       
    %counts mat
    % get rid of row+col 1
    % because dipy records the connection to '0'
    dataRaw(idx).lensMat(1,:) = [] ;
    dataRaw(idx).lensMat(:,1) = [] ;
        
    %% read in vols
    %%%%%%%%%%%%%%%

    tempPotFile = dir(char(strcat(volMatsDir,'*',ids(idx),'*.csv'))) ;    
    % put into struct
    fileName = char(strcat(volMatsDir,tempPotFile.name));
    % read in file
    dataRaw(idx).volsMat = csvread(fileName) ; 
    
    %% read in coordinates
    %%%%%%%%%%%%%%%%%%%%%%
    
    % normalized coor file
    tempPotFile = dir(char(strcat(coorMatsDir,'*',ids(idx),'*.csv'))) ;
    % put into struct
    fileName = char(strcat(coorMatsDir,tempPotFile.name));
    dataRaw(idx).coorMM = csvread(fileName) ; 
    
end

%% fix up some of the NKI data
% in this loop we gonna clean some stuff up
for idx=1:length(dataRaw)
   
    disp(idx)
    
    %% make countVolNorm 
    %%%%%%%%%%%%%%%%%%%%
    
    dataRaw(idx).countVolNormMat = dataRaw(idx).countMat(:,:) ; 
    
    %lets make the streamlinmes / vol mat
    for idx_a=1:length(dataRaw(idx).volsMat)  
        for idx_b=1:length(dataRaw(idx).volsMat)
            
            dataRaw(idx).countVolNormMat(idx_a,idx_b) = ...
                dataRaw(idx).countVolNormMat(idx_a,idx_b) ./ ...
                ( sqrt(dataRaw(idx).volsMat(idx_a,3)) * sqrt(dataRaw(idx).volsMat(idx_b,3)) );

        end
    end
    
    % get rid of nana
    dataRaw(idx).countVolNormMat(isnan(dataRaw(idx).countVolNormMat)) = 0 ;

    %% other stuff
    %%%%%%%%%%%%%%
    
    % make a distance matrix
    dataRaw(idx).distCoorMM = squareform(pdist(dataRaw(idx).coorMM(:,3:5))) ;   
     
    % make vol + len nor
    dataRaw(idx).countVolLenNormMat = dataRaw(idx).countVolNormMat ./ ...
        dataRaw(idx).lensMat ;
    dataRaw(idx).countVolLenNormMat(isnan(dataRaw(idx).countVolLenNormMat)) = 0 ;
    
end

%% add ICV

landR = readtable(strcat(rawDataDir,'/LandRvolumes.csv')) ;
landR.SubjID = strrep(landR.SubjID,'sub-','') ;
landR.Properties.RowNames = landR.SubjID ;
landRsubset = landR(:,{'ICV' 'SupraTentorial' 'SupraTentorialNotVent'...
    'BrainSeg' 'BrainSegNotVent'}) ;

demographicsRaw = join(demographicsRaw , landRsubset, 'Keys', 'RowNames') ;

%% package it all up
% to return from function 

datasetStruct = struct() ;
datasetStruct.dataRaw = dataRaw ;
datasetStruct.demoRaw = demographicsRaw ;

% and save it for future use
save(strcat(rawDataDir,parcellation,'datasetStruct.mat'),'datasetStruct')

