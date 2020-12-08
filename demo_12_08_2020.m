%% demo on Dec 8,2020

% In this demo, we use three photon-efficient algorithms (Shin et al. 2015, Rapp et al. 2017, and Li et al. 2020)
% to process our long-range single-photon lidar data 

clc; clear all; close all;

%% load the data

fprintf('Loading the data... \n \n');
load('./data/tower.mat');

%% test with Shin et al. 2015 and Rapp et al. 2017

addpath('./fcns1');

fprintf('* Reconstruction with Shin et al. 2015...\n\n');
test_Shin;

fprintf('* Reconstruction with Rapp et al. 2017...\n\n');
test_Rapp;

rmpath('./fcns1');

%% test with Li et al. 2020

addpath('./fcns2');

fprintf('* Reconstruction with Li et al. 2020...\n\n');
test_Li;

rmpath('./fcns2');
    
%% plot the reconstruction results    

fprintf('Plot the reconstruction results...\n\n');
figure;
subplot(1,3,1); imagesc(Shin_Z,[0,12]); axis image; axis off; colorbar; colormap jet; title('Shin et al. 2015');
subplot(1,3,2); imagesc(Rapp_Z,[0,12]); axis image; axis off; colorbar; colormap jet; title('Rapp et al. 2017');
subplot(1,3,3); imagesc(Li_Z,[0,12]); axis image; axis off; colorbar; colormap jet; title('Li et al. 2020');
disp('Complete!');






