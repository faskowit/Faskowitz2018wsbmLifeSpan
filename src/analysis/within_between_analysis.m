%% lets analyze within-block vs between-block weighted degree

nSubj = size(fitWSBMAllStruct,2) ;

% with the ya_templateModel
[~,ca_wsbm] = community_assign(templateModel);
nBlocks = size(templateModel.Para.mu,1);

%communities
caIdx = ~~dummyvar(ca_wsbm);

within_strength = zeros([ nBlocks nSubj ]);
between_strength = zeros([ nBlocks nSubj ]);

% iterate over subjects
for idx=1:nSubj
   
    %raw data is already thresholded by the init_thr
    tmpSubjMat = fitWSBMAllStruct(idx).Raw_Data ;
    % get rid of the Nans
    tmpSubjMat(isnan(tmpSubjMat)) = 0 ;
    
    %iterate over blocks/communities
    for jdx=1:nBlocks

        within_strength(jdx,idx) = sum(sum(tmpSubjMat(caIdx(:,jdx),caIdx(:,jdx)))) ;
        between_strength(jdx,idx) = sum(sum(tmpSubjMat(caIdx(:,jdx),~caIdx(:,jdx)))) ;

    end
end