function [blur]=blurcal_3D(samplesize,bin,fovx,fovy,timescale,sigma,width)
blur=zeros(fovx,fovy,timescale);
for a=1:fovx
    for b=1:fovy
        for c=1:timescale
        thetain=acos(cos((a-(fovx+1)/2)*samplesize)*cos((b-(fovy+1)/2)*samplesize));
        blur(a,b,c)=exp(-(thetain/sigma)^2/2-((c-(timescale+1)/2)/width*bin)^2/2);
        end
    end
end
end