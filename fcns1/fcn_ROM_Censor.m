function [T_rom_SM, T_rom_ML, numCensDetect ] = fcn_ROM_Censor(totDetect, I_MAP, Tp, bgndRate, s_fit, romOutput)
% Perform hypothesis testing

[Lr,Lc] = size(totDetect);

T_rom_SM = zeros(Lr,Lc);
T_rom_ML = zeros(Lr,Lc);
numCensDetect = zeros(Lr,Lc);

% Compute peak of pulse shape for ML depth computation
[~,fitPeak] = max(s_fit);

parfor ii=1:Lr
% for ii=1:Lr
    for jj=1:Lc

        data_vec = totDetect{ii,jj};
        if( (romOutput(ii,jj) ~= 0) && ~isempty(data_vec) )
            % when rom has meaningful data
            censDetect = data_vec(  abs(romOutput(ii,jj)-data_vec) < 2*Tp*bgndRate./(I_MAP(ii,jj)+bgndRate));
            numCensDetect(ii,jj) = length(censDetect);
            
            if ~isempty(censDetect)
                T_rom_SM(ii,jj) = mean(censDetect);
                
                % ML Depth from matched filter
                max_temp = max(censDetect);
                min_temp = min(censDetect);
                sub_min = min_temp-1;

                timeBinSig = zeros(max_temp-sub_min,1);
                bins = unique(censDetect);
                timeBinSig(bins-sub_min) = histc(censDetect,bins);

                [xcorr_out, lags] = xcorr(timeBinSig,s_fit);

                % Find height of peak correlation
                mXcor = max(xcorr_out);

                % If multiple lags have same peak, choose one at random
                maxIdces = find(xcorr_out==mXcor);
                mxIndx = maxIdces(randperm(length(maxIdces),1)); 
                
                T_rom_ML(ii,jj) = lags(mxIndx)+fitPeak+sub_min;
                
            end           
        end        
    end
end

end
