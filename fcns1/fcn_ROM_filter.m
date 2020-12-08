function [ romOutput ] = fcn_ROM_filter(totDetect, romSize)
%Compute ROM filter output, use to censor noisy detections for convex depth
%estimation
[Lr,Lc] = size(totDetect);

cell_large = [cell(Lr,romSize) totDetect cell(Lr,romSize)];
cell_large = [cell(romSize,Lc+2*romSize) ; cell_large ; cell(romSize,Lc+2*romSize)];

romOutput = zeros(Lr,Lc);


parfor ii = 1:Lr
%for ii = 1:Lr
    for jj=1:Lc
        dats = cell_large((ii:ii+2*romSize),(jj:jj+2*romSize));
        cntrPx = 2*romSize^2+2*romSize+1;
        dats_vector = [dats(1:cntrPx-1),dats(cntrPx+1:end)]';   
        raw_time = cell2mat(dats_vector);%(~cellfun('isempty',dats_vector)));        
        filt_val = median(raw_time);
        romOutput(ii,jj) = filt_val;
    end
end
romOutput(isnan(romOutput)) = 0;


end

