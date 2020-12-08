%% A deconvolution Reconstruction Scheme
% A 3D deconvolution scheme modified from SPIRAL-TAP by Z.H.Harmany et al.

%% Preprocess the uncensored data and Set some algorithm parameters

% translate the type of the data from cell to 3D matrix
totDetect = data_processed;                                      % read the data
[lr,lc] = size(totDetect);                                       % size of the picture
bin=1*10^(-9);                                                   % temporal resolution of the 3D matrix is 1ns, 
                                                                 % unit: s
c_v = 3e8;

y = preprocess_cell2matrix(MM, mm, lr, lc, bin, totDetect);      % turn the cell totDetect to a 3D matrix y 
                                                                 % with 1ns bin resolution

% ganerate the imaging forward model 
% according to the calibration of the imaging setup 
width = 10^(-9);                                                % the whole system jitter is 1ns
sigma = 1.78e-5;                                                % the size of the light spot 
samplesize = 1.12e-5;                                           % the size of the FoV
                                                                 % over 1400m range, 
                                                                 % the radii of the two above are 0.025m and 0.015m 
                                                                 
blur=blurcal_3D(samplesize,bin,3,3,21,sigma,width);              % generate the 3D spatiotemporal kernel       
A = @(x) convn(x,blur,'same'); AT = @(x) convn(x,blur,'same');   % construct the imaging forward model                              
        
% set parameters :
tautv = 0.1;
miniter = 1;
maxiter = 2;
stopcriterion = 3;
tolerance = 1e-8;
verbose = 10;        

    

%% Reconstruction with a convex solver modified from SPIRALTAP

[recoveredSPIRALtv, iterationsSPIRALtv, objectiveSPIRALtv,...
reconerrorSPIRALtv, cputimeSPIRALtv] ...
         = SPIRALTAP3D_v1(y,A,tautv,...
            'noisetype','gaussian',...
            'penalty','TV',...
            'AT',AT,...    
            'maxiter',maxiter,...
            'Initialization',y,...
            'miniter',miniter,...
            'stopcriterion',stopcriterion,...
            'tolerance',tolerance,...
            'monotone',1,...
            'saveobjective',1,...
            'savereconerror',1,...
            'savecputime',1,...
            'savesolutionpath',0,...
            'truth',y,...
            'verbose',verbose);

%% filter the depth map with the reflective map

[Li_Alpha,Li_Z] = max(recoveredSPIRALtv,[],3);
Li_Z(Li_Alpha<0.1*sum(Li_Alpha(:))/lr/lc) = 0;
Li_Z = Li_Z*c_v*bin/2;

