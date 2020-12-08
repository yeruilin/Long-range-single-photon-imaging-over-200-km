function [ LIDAR_Alpha, LIDAR_Z ] = fcn_ML(SBR,spad)
% Compute conventional LIDAR ML estimates:
%   reflectivity: normalized photon count
%   depth: delay satisfying log-matched filter

c = 3e8;
bin_size = 26e-12; % the temporal res

pulse_bins = (440e-12) / bin_size;
pulse = normpdf(1:8*pulse_bins,(8*pulse_bins-1)/2,pulse_bins/2);


totDetect=cell(size(spad,1),size(spad,2));

for i=1:size(spad,1)
    for j=1:size(spad,2)
        tempvect=[];
        for k=1:size(spad,3)
            photon_counts=round(spad(i,j,k));
            binvect=k*ones(photon_counts,1);
            tempvect=[tempvect;binvect];
        end
        totDetect{i,j}=sort(tempvect);
    end
end

numTotDetect = cellfun('length',totDetect);
numFrames = max(numTotDetect(:));

LIDAR_Alpha = numTotDetect/numFrames;
LIDAR_Z = 0.5*c*fcn_match_filt(totDetect, pulse)*bin_size;


end