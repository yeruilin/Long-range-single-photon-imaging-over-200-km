function [ Ncl_rule ] = fcn_Ncl_rule2( numCl,spxlSizes,bgndRate,numFrames,Tr,Twind,tau_FA)

maxSuperpixelSize = max(spxlSizes);
maxPoissLambda = ceil(bgndRate*numFrames*(2*maxSuperpixelSize+1)^2)+1;
maxNumDetect = ceil(maxPoissLambda+10*sqrt(maxPoissLambda));

Twind_scaled = Twind/Tr;

noiseLambda = 1:maxPoissLambda;

Ncls = 1:numCl;
noiseDets = 1:maxNumDetect;

%% Compute Probability for each Ncl
[ns, qs] = meshgrid(noiseDets,Ncls);

Pr_oneclust = cdf('beta',Twind_scaled,qs,ns-qs+1);
Pr_oneclust(isnan(Pr_oneclust)) = 0; 
Pr_cond = 1-(1-Pr_oneclust).^(ns-qs);

Pr_uncond = zeros(numCl,maxPoissLambda);

parfor ii = Ncls
    for jj = 1:maxPoissLambda
        % only need to calculate pr(poisson) for limited range of
        % possibilities
        lambda_n = noiseLambda(jj);
        rangePoiss = max(floor(lambda_n-10*sqrt(lambda_n)),1):ceil(lambda_n+10*sqrt(lambda_n));
        
        Pr_poiss = poisspdf(rangePoiss,lambda_n);
        Pr_uncond(ii,jj) = Pr_poiss*Pr_cond(ii,rangePoiss)';%
    end
end

%% Choose Threshold Tau_FA
Pr_NoiseClust = Pr_uncond';
idx = zeros(numCl,1);
Ncl_rule = 2*ones(maxPoissLambda,1);
Pr_Ncl_noise = Pr_NoiseClust(:,1);

for ii = 1:numCl
    if size(find(Pr_NoiseClust(:,ii)<tau_FA,1,'last'),1)<1
        idx(ii)=1;
        disp('something strange\n');
    else
        idx(ii) = find(Pr_NoiseClust(:,ii)<tau_FA,1,'last');
    end                                         %HX
    %idx(ii) = find(Pr_NoiseClust(:,ii)<tau_FA,1,'last');
    Ncl_rule(idx(ii)+1:end,1) = (ii+2)*ones(maxPoissLambda-idx(ii),1);

    switch ii
        case numCl
            Pr_Ncl_noise(idx(ii)+1:end,1) = 0.1*ones(maxPoissLambda-idx(ii),1);
        otherwise
            Pr_Ncl_noise(idx(ii)+1:end,1) = Pr_NoiseClust(idx(ii)+1:end,ii+1);
    end
end

end