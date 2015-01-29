%% startup
clear all, close all, clc

path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/figures';
cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private');

%% linker proberties
linker_mean = [8.18 10.16 10.65 11.13 14.53]/10 ;   % nm, BMOE, BMH mPDM, pPDM, BM
linker_sigma = [0.75 2.41 0.55 0.52 1.51] / 10; % nm, fluctuations of maleimide linker
linker_names = {'BMOE', 'BMH', 'mPDM', 'pPDM', 'BM'};

%% 
%data = load('/Users/jonasfunke/Documents/FRET_STAGE/Enzyme/2014-10-14_linker-comparison/Linker-comparison.mat');
data = load('/Users/jonasfunke/Dropbox/POSITIONING/2014-10-14_linker-comparison/data_02.mat');
%%
ke = 0.17; %/uM/h
c_linker = 50; % uM

tic

model_BMH = fittype( @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(2), linker_sigma(2) ) , ...
    'independent', {'d_mean'});
fit_ssLeash = fit(data.experiment.d', data.s_control, model_BMH, 'startpoint', [1.97], 'Lower', [0 ], 'Upper', [2 ]) 

toc
%%
xplot_ssLeash = [0:0.5:7]';
yplot_ssLeash = fit_ssLeash(xplot_ssLeash);

%%

model_BMOE = fittype( @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(1), linker_sigma(1) ) , ...
    'independent', {'d_mean'});
fit_BMOE = fit(data.d16, data.yield16(:,1), model_BMOE, 'startpoint', [0.619], 'Lower', [0 ], 'Upper', [2 ]) 

model_BMH = fittype( @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(2), linker_sigma(2) ) , ...
    'independent', {'d_mean'});
fit_BMH = fit(data.d16, data.yield16(:,2), model_BMH, 'startpoint', [0.619], 'Lower', [0 ], 'Upper', [2 ]) 

model_pPDM = fittype( @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(4), linker_sigma(4) ) , ...
    'independent', {'d_mean'});
fit_pPDM = fit(data.d16, data.yield16(:,3), model_pPDM, 'startpoint', [0.619], 'Lower', [0 ], 'Upper', [2 ]) 

model_BM = fittype( @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(5), linker_sigma(5) ) , ...
    'independent', {'d_mean'});
fit_BM = fit(data.d16, data.yield16(:,4), model_BM, 'startpoint', [0.619], 'Lower', [0 ], 'Upper', [2 ]) 
toc


%%
fit_BMOE.d_sigma
fit_BMH.d_sigma
fit_pPDM.d_sigma
fit_BM.d_sigma
fit_ssLeash.d_sigma

%%
yplot1 = fit_BMOE(data.d16);
yplot2 = fit_BMH(data.d16);
yplot3 = fit_pPDM(data.d16);
yplot4 = fit_BM(data.d16);
%%
cc = [ 0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

close all
cur_fig = figure();

plot(  data.d16, yplot1, '-', 'Color', cc(1,:)), hold on
plot( data.d16, data.yield16(:,1), '.', 'Color', cc(1,:))
plot(  data.d16, yplot2, '-', 'Color', cc(2,:)), hold on
plot( data.d16, data.yield16(:,2), '.', 'Color', cc(2,:))
plot(  data.d16, yplot3, '-', 'Color', cc(3,:)), hold on
plot( data.d16, data.yield16(:,3), '.', 'Color', cc(3,:))
plot(  data.d16, yplot4, '-', 'Color', cc(4,:)), hold on
plot( data.d16, data.yield16(:,4), '.', 'Color', cc(4,:))

plot( xplot_ssLeash, yplot_ssLeash, '-', 'Color', cc(5,:))
plot( data.experiment.d', data.s_control, '.', 'Color', cc(5,:))

%%

save('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private/fit_thiol_fluctuations_data.mat')

%% 
ke = 0.17; %/uM/h
c_linker = 50; % uM


model_global = fittype( @(d_sigma, d_mean, linker_mean, linker_sigma) compute_normalized_yield_global(ke, c_linker, d_mean, d_sigma, linker_mean, linker_sigma ) , ...
    'independent', {'d_mean', 'linker_mean', 'linker_sigma'});




