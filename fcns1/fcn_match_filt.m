function [ T_match ] = fcn_match_filt( totDetect,  s_fit)
%Apply matched filter for depth estimation

fprintf('Start: Matched Filtering \n');

[Lr,Lc] = size(totDetect);

T_match = zeros(Lr,Lc);

% Compute peak of pulse shape for ML depth computation
[~,fitPeak] = max(s_fit);

parfor jj=1:Lc
    for ii=1:Lr
        data_vec = totDetect{ii,jj};
            
        if ~isempty(data_vec)
            % ML Depth from matched filter
            max_temp = max(data_vec);
            min_temp = min(data_vec);
            sub_min = min_temp-1;

            timeBinSig = zeros(max_temp-sub_min,1);
            bins = unique(data_vec);
            timeBinSig(bins-sub_min) = histc(data_vec,bins);

            [xcorr_out, lags] = xcorr(timeBinSig,s_fit);

            % Find height of peak correlation
            mXcor = max(xcorr_out);

            % If multiple lags have same peak, choose one at random
            maxIdces = find(xcorr_out==mXcor);
            mxIndx = maxIdces(randperm(length(maxIdces),1)); 

            T_match(ii,jj) = lags(mxIndx)+fitPeak+sub_min;

        end           

    end
end
fprintf('End: Matched Filtering \n \n');

end

