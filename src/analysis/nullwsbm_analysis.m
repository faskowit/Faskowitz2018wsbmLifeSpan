
nSubj = size(fitWSBMAllStruct,2) ;

%% calculate how each observed block matrix deviates from model block matrix
% the null wsbm is different than the emperical 
% the two nulls below are the same...because we are using harsh_mu...
% so there is actually no inference to be done

% get the harsh_mu, with binary probabilities
[ ~ , harsh_mu ] = make_WSBM_prior(templateModel , 1) ;

% make NULL model input
% (need to have at least one trial...but it wont do anything) 
nullModelStruct = indivModelInputStruct ;
nullModelStruct.mu_0 = harsh_mu ;
nullModelStruct.numTrials = 1 ;
nullModelStruct.muMaxIter = 0 ;
nullModelStruct.mainMaxIter = 0 ;
nullModelStruct.verbosity = 0 ;
% make this struct into a cell list that the wsbm likes
a = struct2nv(nullModelStruct) ;
b = struct2cell(nullModelStruct) ;
c = [ a b ]';
nullModelInput = c(:)' ;

fitWSBMAllnull = struct() ;

for idx=1:nSubja
   
    disp(idx)
    
    [~,fitWSBMAllnull(idx).nullModel] = wsbm(fitWSBMAllStruct(idx).Raw_Data,...
        fitWSBMAllStruct(idx).Model(1).R_Struct.R,...
        nullModelInput{:}) ;
 
end

% % lets look at this real quick
% tmpVec = zeros([ nSubj 3 ]);
% 
% for idx=1:nSubj
% 
%     tmpVec(idx,1) = tmp(idx).nullModel.Para.LogEvidence ;
%     tmpVec(idx,2) = fitWSBMAllStruct(idx).Model(5).Para.LogEvidence ;
%     tmpVec(idx,3) = datasetDemo_uncond{idx,'age'} ;
%     
% end

% [~,tmpModel] = wsbm(fitWSBMAllStruct(1).Raw_Data,fitWSBMAllStruct(1).Model(5).R_Struct,...
%     modelInputs{:}) ;

%% calculate proportion of group connectivity...

nBlocks = size(templateModel.Para.mu,1);

%communities
ca = community_assign(templateModel.Para.mu);
caIdx = ~~dummyvar(ca(:,2));

observedE_Pct = zeros([nBlocks nBlocks nSubj]); 

% iterate over subjects
for idx=1:nSubj
   
    tmpSubjMat = fitWSBMAllStruct(idx).Raw_Data ;
    % get rid of the Nans
    tmpSubjMat(isnan(tmpSubjMat)) = 0 ;
    
    %iterate over blocks/communities
    for jdx=1:nBlocks
        
%         [a,b] = size(tmpSubjMat(caIdx(:,jdx),:));
%         size_totCom = a * b ;
%         
%         e_totCom = sum(sum(tmpSubjMat(caIdx(:,jdx),:) > 0));
%         w_totCom = sum(sum(tmpSubjMat(caIdx(:,jdx),:))) ;
        
        for kdx=1:nBlocks
        
            if kdx < jdx
                continue
            end
            
            [c,d] = size(tmpSubjMat(caIdx(:,jdx),caIdx(:,kdx)));
            size_com = c * d ;
            
            %w_com = tmpSubjMat(caIdx(:,jdx),caIdx(:,kdx)) ;
            e_com = sum(sum(tmpSubjMat(caIdx(:,jdx),caIdx(:,kdx)) > 0)) ;
            
            e_com_pct = e_com / size_com ;
            
            observedE_Pct(jdx,kdx,idx) = e_com_pct ; 
            
        end
    end
end

%% correlate observed with wsbm 

obsCormodel_Vec_1 = zeros([ nSubj 1]);
obsCormodel_Vec_2 = zeros([ nSubj 1]);

ind = ~~triu(ones(nBlocks));

for idx=1:nSubj
    
    tmpSubjMat = observedE_Pct(:,:,idx) ;
    obsCormodel_Vec_1(idx) = corr(fitWSBMAllStruct(idx).Model(5).Para.predict_e,...
        tmpSubjMat(ind));
    
    obsCormodel_Vec_2(idx) = corr(fitWSBMAllStruct(idx).Model(5).Para.predict_e,...
        fitWSBMAllnull(idx).nullModel.Para.predict_e);  
       
end







