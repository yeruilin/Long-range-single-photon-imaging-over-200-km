% Unmixing Processing      
% 1) Windowing
% 2) Regularized Refl. Estimation
% 3) Superpixel Formation
% 4) Windowing of superpixels
% 5) Regularized Reflectivity Re-estimation
% 6) ROM Censoring 
% 7) Regularized Depth estimation

%% Fixed Lidar system parameters
c_v = 3e8;
bin_size = 1e-12; % the temporal res
pulse_bins = (1e-9) / bin_size;
pulse = normpdf(1:8*pulse_bins,(8*pulse_bins-1)/2,pulse_bins/2);


%% Read and preprocess the data

totDetect = data_processed;
mean_signal_photons = signal_per_pixel;
SBR = sbr;

numTotDetect = cellfun('length',totDetect);
numFrames = pulse_num;
NrB = mean_signal_photons/SBR;                    
bgndRate = NrB/numFrames;
Twind = 2*pulse_bins;
Tr = MM-mm;
[Lr,Lc] = size(numTotDetect);

parfor iii = 1:Lr
    for jjj = 1:Lc
        totDetect{iii,jjj} = sort(totDetect{iii,jjj}-mm)';
    end
end


%% Fixed Algorithm Parameters

%tv
beta_alpha = 0.1; %1
beta_z = 0.1;

% Superpixels
nmRng = 0.05;

% Alpha Regularization
max_iter_alpha = 40;
pen_type = 'tv'; 
AT = @(x) x; A = @(x) x;

% ROM Censoring
romSize = 2;    
nearest_k = 7;  

% Z Regularization
max_iter_z = 20; 
ainit_z = 0.05; 

tau_FA = 0.01;          % Probability of false acceptance threshold  0.01 
maxDsp = 1;             % Max superpixel size (determines number of iterations)


%% Reconstruction

mainSuffIndex = zeros(Lr,Lc);
numCombDetect = zeros(Lr,Lc);
combDetect = cell(Lr,Lc);

Ncl_rule = fcn_Ncl_rule2( 80, 0:maxDsp,bgndRate,numFrames,Tr,Twind,tau_FA);%1 to 0

noise_window = bgndRate*(Twind/Tr);

for spxlSize = 0:maxDsp
    % Step 1: Form Superpixels
    switch spxlSize
        case 0  % Single Pixel case - data is raw data
            superDetect = totDetect;
            numSuperpixel = ones(Lr,Lc);
            %beta_alpha = Beta_Alphas(3);
            
        otherwise   % Form superpixels
            [ superDetect, numSuperpixel ] = ...
                fcn_superpixel(totDetect, numTotDetect, bestUnmixingAlpha, spxlSize, nmRng, mainSuffIndex );
            %beta_alpha = Beta_Alphas(5);
    end
    
    % Step 2: Window Superpixels
    [ ~, numWindDetect, windDetect ] = fcn_windowing(superDetect,Twind);

    % Step 3: Combine data
    Ncl = Ncl_rule(ceil(NrB*numSuperpixel)+1);
    mainSuffIndex(numWindDetect>=Ncl)=1;
    
    tempSuffIndex = zeros(Lr,Lc);
    tempSuffIndex(numWindDetect>=Ncl)=1;
    
    numCombDetect(tempSuffIndex==1) = numWindDetect(tempSuffIndex==1)./numSuperpixel(tempSuffIndex==1);
    combDetect(tempSuffIndex==1) = windDetect(tempSuffIndex==1); 
    
    % Step 4: Regularize reflectivity
    insuffNumInit = zeros(Lr,Lc);
    insuffNumInit(mainSuffIndex==0) = numWindDetect(mainSuffIndex==0)./numSuperpixel(mainSuffIndex==0);
    numInit = numCombDetect + insuffNumInit;
    
    reflData = numInit;

    
    scale_factor = max(numSuperpixel(:));
    
    C_tv = SPIRALTAP(round(scale_factor*reflData),A,beta_alpha, ...
        'Initialization',scale_factor*reflData,...
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
    C_tv = C_tv/scale_factor;

    bestUnmixingAlpha = max(C_tv/numFrames-noise_window,0);
    
    
end
for ii = 1:Lr
    for jj = 1:Lc
        dep_comD(ii,jj) = mean(combDetect{ii,jj});
    end
end
dep_comD(isnan(dep_comD)) = 0;


%%

romSize = 2; 
nearest_k = 7; 

% Step 5: Additiomal ROM Censoring
romOutput = fcn_ROM_filter(combDetect, romSize);
[~, T_ML_Unmixing, ~ ] = fcn_ROM_Censor(combDetect, bestUnmixingAlpha, pulse_bins, bgndRate, pulse, romOutput);

% Step 6: Depth - Initializes with kNN to speed up regularization
T_Data = T_ML_Unmixing;
T_Init = fcn_kNN_Init( bestUnmixingAlpha, T_ML_Unmixing, nearest_k );

T_tv_Unmixing = SPIRALTAP(dep_comD*c_v/2*bin_size,A,beta_z, ...
    'noisetype','gaussian', ...
    'penalty',pen_type, ...
    'maxiter',max_iter_z,...
    'Initialization',dep_comD*c_v/2*bin_size,... 
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

bestUnmixingZ = T_tv_Unmixing;

Rapp_Alpha = bestUnmixingAlpha;
%Rapp_Z = bestUnmixingZ;
Rapp_Z = dep_comD*c_v/2*bin_size;

