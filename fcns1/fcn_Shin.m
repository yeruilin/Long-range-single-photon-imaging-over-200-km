function [ bestShinAlpha, bestShinZ ] = fcn_Shin(mean_signal_photons, SBR, spad, beta_z_in)
% Shin TCI Processing
% Reflectivity estimation assuming detection counts are binomial with a
% subtracted background rate. 
% ROM censoring to attempt to remove noise detections. 
% Depth estimation from regularized log-matched filter.

% Shin_Alphas = zeros(Lr,Lc,numBeta);
% MSE_TV_Shin_Alphas = zeros(numBeta,1);
% 
% Shin_Zs = zeros(Lr,Lc,numBeta);
% RMSE_TV_Shin_Zs = zeros(numBeta,1);

%% Fixed Algorithm Parameters

%tv
beta_alpha=0.12;
beta_z=beta_z_in;

% Alpha Regularization
max_iter_alpha = 40;
pen_type = 'tv'; 
AT = @(x) x; A = @(x) x;

% ROM Censoring
romSize = 2;

% Z Regularization
max_iter_z = 20;%100; 
ainit_z = 0.05; 

%%
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
NrB = mean_signal_photons/SBR;                    %total background counts per pixel
bgndRate = NrB/numFrames;


% Step 1: Reflectivity
reflData = numTotDetect;
% I_ML_Shin = max(log(numFrames./(numFrames-reflData))-bgndRate,0);
% MSE_ML_Shin_Alpha =  10*log10(mean((Alpha_true(:)-I_ML_Shin(:)).^2));

fprintf('Start: SPIRAL-Alpha \n');
% parfor balpha = 1:numBeta

    C_tv = SPIRALTAP(reflData,A,beta_alpha, ...
        'noisetype','poisson', ... 
        'penalty','tv', ...
        'maxiter',max_iter_alpha,... 
        'AT',AT,...
        'monotone',1,...
        'stopcriterion',3,...
        'tolerance',1e-12,...
        'saveobjective',1,...
        'savereconerror',0,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'alphainit',0.2,...
        'alphaaccept',1e80,...
        'verbose',0);    
    
    bestShinAlpha = max(log(numFrames./(numFrames-C_tv))-bgndRate,0);
%     Shin_Alphas(:,:,balpha) = I_TV_Shin;
%     MSE_TV_Shin_Alphas(balpha) = 10*log10(mean((Alpha_true(:)-I_TV_Shin(:)).^2));
% end
fprintf('End: SPIRAL-Alpha \n \n');

% bestShinBetaAlpha = balpha;
% bestShinAlpha = Shin_Alphas(:,:,bestShinBetaAlpha);

% Step 2: ROM Censoring
romOutput = fcn_ROM_filter(totDetect, romSize);
[~, T_ML_Shin, ~ ] = fcn_ROM_Censor(totDetect, bestShinAlpha, pulse_bins, bgndRate, pulse, romOutput);

% Step 3: Depth
T_Data = T_ML_Shin;
T_Init = romOutput;

% Z_ML_Shin = T_ML_Shin*ttd;
% RMSE_ML_Shin_Z =  sqrt(mean((Z_true(:)-Z_ML_Shin(:)).^2));

fprintf('Start: SPIRAL-Z \n');
% parfor bz = 1:numBeta
% bz = 12;
    %beta_z = Beta_Zs(bz);
    T_tv_Shin = SPIRALTAP(T_Data,A,beta_z, ...
        'noisetype','gaussian', ...
        'penalty',pen_type, ...
        'maxiter',max_iter_z,...
        'Initialization',T_Init,... 
        'monotone',1,...
        'miniter',1,...
        'AT',AT,...
        'stopcriterion',4,...
        'tolerance',1e-6,...
        'alphainit',ainit_z,...
        'alphaaccept',1e80,...
        'logepsilon',1e-10,...
        'saveobjective',1,...
        'savereconerror',1,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'truth',T_Init,...
        'reconerrortype',0,...
        'verbose',0);

    bestShinZ = T_tv_Shin*c/2*bin_size;
%     Shin_Zs(:,:,bz) = Z_tv_Shin;
%     RMSE_TV_Shin_Zs(bz) = sqrt(mean((Z_true(:)-Z_tv_Shin(:)).^2));
% end
fprintf('End: SPIRAL-Z \n \n');

% bestShinBetaZ = bz;
% bestShinZ = Shin_Zs(:,:,bestShinBetaZ);

disp('Shin TCI completed');

end

