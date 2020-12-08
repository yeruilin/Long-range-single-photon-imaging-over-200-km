function [ superDetect, numSuperpixel ] = fcn_superpixel(totDetect, numTotDetect, I_tv, sprPxlSize, nmRng, suffIndex )
%FCN_SUPERPIXEL creates expanded detection vectors by borrowing detections
%   from similar neighboring pixels if the original detection vector had an
%   insufficient number of signal detections.
%   Similar neighbors are chosen from a set of pixels of distance at most
%   sprPixelSize pixels away from the center and with reflectivity estimate
%   at most nmRng from the center's reflectivity.

[Lr,Lc] = size(I_tv);

pcTile = prctile(I_tv(:),[2,98]);
scale_image = (I_tv-pcTile(1))/(pcTile(2)-pcTile(1));
scale_image(scale_image > 1) = 1;
scale_image(scale_image < 0) = 0;

% Parameters
preDetect = cell(Lr+2*sprPxlSize,Lc+2*sprPxlSize);
preDetect(1+sprPxlSize:Lr+sprPxlSize,1+sprPxlSize:Lc+sprPxlSize) = totDetect;
superDetect = cell(Lr,Lc);%
numSuperpixel = ones(Lr,Lc);

pnumDetect = -10*ones(Lr+2*sprPxlSize,Lc+2*sprPxlSize);
pnumDetect(1+sprPxlSize:Lr+sprPxlSize,1+sprPxlSize:Lc+sprPxlSize) = scale_image;


parfor jj = 1:Lc
%for jj = 1:Lc
    if ~mod(jj,floor(Lr*Lc/10))
        fprintf('*');
    end
    for ii = 1:Lr
        if suffIndex(ii,jj)==0
            pxRefl = pnumDetect(ii+sprPxlSize,jj+sprPxlSize);

            dtctWndw = preDetect(ii:ii+2*sprPxlSize,jj:jj+2*sprPxlSize);
            numWndw = pnumDetect(ii:ii+2*sprPxlSize,jj:jj+2*sprPxlSize);
            goodRefl = abs(pxRefl-numWndw) <= nmRng;
            goodCell = cell(2*sprPxlSize+1);
            goodCell(goodRefl==1) = dtctWndw(goodRefl==1);
            superDetect{ii,jj} = sort(cell2mat(goodCell(:)));
            numSuperpixel(ii,jj) = sum(goodRefl(:)==1);
        end
    end
end

end