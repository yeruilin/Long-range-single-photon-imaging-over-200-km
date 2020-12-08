function [mlWindDepth, numWindDetect, windDetect] = fcn_windowing(totalDetect,Twind)
%FCN_FIXED_WINDOW uses overlapping windows to find the best window for a
%given window size
[Lr,Lc]= size(totalDetect);

mlWindDepth = zeros(Lr,Lc); % Distance of max correlation
numWindDetect = zeros(Lr,Lc); % Number of valid detections
windDetect = cell(Lr,Lc);

Ltvect = Lr*Lc;

parfor kk = 1:Ltvect   
% for kk = 1:Ltvect      
    receiveTT = totalDetect{kk};
    numDets = numel(receiveTT);
       
    if numDets > 1 % If pixel has at least 1 detection
        % Compute time difference between detections
        detDiff = diff(receiveTT);
        numWind = zeros(numDets,1);
        
        % For each det, check whether subsequent detections occur within
        % window duration
        for ii = 1:numDets
            cumDiff = cumsum(detDiff(ii:end));
            numWind(ii) = sum(cumDiff<Twind);
        end
        
%         % For each det, check whether subsequent detections occur within
%         % window duration        
%         cumDiffMat = zeros(numDets);
%         cumDiffMat(:,2:end) = repmat(cumsum(detDiff'),numDets,1);
%         cDiffMat = cumDiffMat-cumDiffMat'-tril(ones(numDets),0);
%         cDiffMat(cDiffMat>Twind)=-1;
%         
%         numWind = sum(cDiffMat>=0,2);
                     
        % Find number of detections in window with max
        mX = max(numWind);
        
        % If multiple windows have same max number, choose one at random
        maxIdces = find(numWind==mX);
        mxIndx = maxIdces(randperm(length(maxIdces),1)); 

        % Return size of cluster (num detections) and position of cluster
        numWindDetect(kk) = 1+mX;
        windDetect{kk} = receiveTT(mxIndx:mxIndx+mX);
        mlWindDepth(kk) = mean(windDetect{kk});
        
    elseif numDets == 1
        numWindDetect(kk) = 1;
        windDetect{kk} = receiveTT;
        mlWindDepth(kk) = receiveTT;
    end
end


end

