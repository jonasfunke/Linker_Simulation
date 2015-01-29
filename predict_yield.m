%% startup
clear all, close all, clc

%path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/figures';
%cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private');

%% linker proberties
linker_mean = [8.18 10.16 10.65 11.13 14.53]/10 ;   % nm, BMOE, BMH mPDM, pPDM, BM
linker_sigma = [0.75 2.41 0.55 0.52 1.51] / 10; % nm, fluctuations of maleimide linker
linker_names = {'BMOE', 'BMH', 'mPDM', 'pPDM', 'BM'};

%%
BMOE.mean = linker_mean(1);
BMOE.std = linker_sigma(1);
BME.name = 'BMOE';

BMH.mean = linker_mean(2);
BMH.std = linker_sigma(2);
BMH.name = 'BMH';

Linker.BMH = BMH;
Linker.BMOE = BMOE;
%% 
%data = load('/Users/jonasfunke/Documents/FRET_STAGE/Enzyme/2014-10-14_linker-comparison/Linker-comparison.mat');
%data = load('/Users/jonasfunke/Dropbox/POSITIONING/2014-10-14_linker-comparison/data_02.mat');
%%
ke = 0.17; %/uM/h
c_linker = 10; % uM

yield_BMOE = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(1), linker_sigma(1) );
yield_BMH = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(2), linker_sigma(2) );
yield_BM = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(5), linker_sigma(5) );

%%
n = [10:50]';
d47 = 0.335.*dsDNA_mean_end_to_end_distances_semiflexible_polymer(n, 4)*10; % in [A]
d15 = (15/46.5).*d47(:,1);

sigma = 0.1:0.1:1;
y1 = zeros(length(d15), length(sigma));
y2 = zeros(length(d15), length(sigma));
y3 = zeros(length(d15), length(sigma));

for i=1:length(sigma)
    y1(:,i) = yield_BMOE(sigma(i), d15/10);
    y2(:,i) = yield_BMH(sigma(i), d15/10);
    y3(:,i) = yield_BM(sigma(i), d15/10);
end



close all
plot(d15, y1, d15, y2, d15, y3)








