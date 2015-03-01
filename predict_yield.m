%% startup
clear all, close all, clc

path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/';
%cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private');

path_out = ''
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

xplot = 0:0.1:60;

sigma = 0.1:0.1:1;
y1 = zeros(length(xplot), length(sigma));
%y2 = zeros(length(d15), length(sigma));
%y3 = zeros(length(d15), length(sigma));

for i=1:length(sigma)
    disp(['sigma = ' num2str(sigma(i)) ' nm'])
    tic
    y1(:,i) = yield_BMOE(sigma(i), xplot/10);
  %  y2(:,i) = yield_BMH(sigma(i), d15/10);
  %  y3(:,i) = yield_BM(sigma(i), d15/10);
    toc
end

%%
save([path_out 'data' filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_data.mat'])

%%
close all


cc = varycolor(length(sigma));

thiol1 = @(x, sigma) pdf('Normal', x, 0, sigma*10);
thiol2 = @(x, sigma) pdf('Normal', x, 30, sigma*10);
xplot_thiol = -30:0.1:80;

ythiol = zeros(length(xplot_thiol), length(sigma)*2);

for i=1:length(sigma)
    subplot(1, 2, 1)
    plot(xplot, y1(:,i), 'Color', cc(i,:)), hold on
    
    subplot(1, 2, 2)
    plot(xplot, thiol1(xplot, sigma(i)), 'Color', cc(i,:)), hold on
    plot(xplot, thiol2(xplot, sigma(i)), 'Color', cc(i,:)), hold on
   % plot(xplot, thiol1(xplot, sigma(i)).*thiol2(xplot, sigma(i)), 'Color', cc(i,:)), hold on

   ythiol(:, 2*(i-1)+1) = [thiol1(xplot_thiol, sigma(i))'];
   ythiol(:, 2*(i-1)+2) = [thiol2(xplot_thiol, sigma(i))'];
    
end

%%
bla = [xplot', y1];
bla_thiol = [xplot_thiol', ythiol];




%%
close all
subplot(3,1,1)
plot(d15, y1)

subplot(3,1,2)
plot(d15, y2)

subplot(3,1,3)
plot(d15, y3)


%%
ke = 0.17; %/uM/h
c_linker = 10; % uM

yield_BMOE = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(1), linker_sigma(1) );

xx = d15(1):0.1:d15(end)
yy = yield_BMOE(sigma(1), xx/10);

%%
BNOE_out_10uM = [d15, y1];

%xx2 = d15(1):0.1:d15(end);

yy2 = zeros(length(xx), 3);

yy2(:,2) = spline(d15, y1(:,5), xx);
yy2(:,3) = spline(d15, y1(:,10), xx);

%%
out_fine = [xx', yy', yy2(:,2:3)];

%%
cc = varycolor(length(sigma));

thiol1 = @(x, sigma) pdf('Normal', x, 0, sigma*10);
thiol2 = @(x, sigma) pdf('Normal', x, 30, sigma*10);
xplot = -30:0.1:80;
close all
for i=1:length(sigma)
    subplot(1, 2, 1)
    plot(d15, y1(:,i), 'Color', cc(i,:)), hold on
    
    subplot(1, 2, 2)
    plot(xplot, thiol1(xplot, sigma(i)), 'Color', cc(i,:)), hold on
    plot(xplot, thiol2(xplot, sigma(i)), 'Color', cc(i,:)), hold on
   % plot(xplot, thiol1(xplot, sigma(i)).*thiol2(xplot, sigma(i)), 'Color', cc(i,:)), hold on

    
end


%%
cc = varycolor(length(sigma));

close all
hold all
for i=1:length(sigma)
    plot(d15, y1(:,i), 'Color', cc(i,:))
    plot(d15, y2(:,i), 'Color', cc(i,:))
    plot(d15, y3(:,i), 'Color', cc(i,:))
end


%% 2 mM Linker concentration
ke = 0.17; %/uM/h
c_linker = 2000; % uM

yield_BMOE = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(1), linker_sigma(1) );
yield_BMH = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(2), linker_sigma(2) );
yield_BM = @(d_sigma, d_mean) compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean(5), linker_sigma(5) );

y1_2mM = zeros(length(d15), length(sigma));
y2_2mM = zeros(length(d15), length(sigma));
y3_2mM = zeros(length(d15), length(sigma));

for i=1:length(sigma)
    disp(['sigma = ' num2str(sigma(i)) ' nm'])
    tic
    y1_2mM(:,i) = yield_BMOE(sigma(i), d15/10);
    y2_2mM(:,i) = yield_BMH(sigma(i), d15/10);
    y3_2mM(:,i) = yield_BM(sigma(i), d15/10);
    toc
end


%%

cur_fig = figure();
hold all
for i=1:length(sigma)
    subplot(3,1,1)
    plot(d15, y1(:,i), 'Color', cc(i,:)), hold on
    plot(d15, y1_2mM(:,i), 'Color', cc(i,:)), hold on
    title('BMOE-Linker')
    
    subplot(3,1,2)
    plot(d15, y2(:,i), 'Color', cc(i,:)), hold on
    plot(d15, y2_2mM(:,i), 'Color', cc(i,:)), hold on
    title('BMH-Linker')
    
    subplot(3,1,3)
    plot(d15, y3(:,i), 'Color', cc(i,:)), hold on
    plot(d15, y3_2mM(:,i), 'Color', cc(i,:)), hold on
    title('BM-Linker')

end

print(cur_fig, '-dtiff', '-r500', [path_out 'figures' filesep 'Prediction_yield.tif'])

%%

close all
plot(d15, y1, d15, y1_2mM)


%%
close all
subplot(3,1,1)
plot(d15, y1_2mM)

subplot(3,1,2)
plot(d15, y2_2mM)

subplot(3,1,3)
plot(d15, y3_2mM)

%% save
mkdir([path_out 'data' filesep])
save([path_out 'data' filesep datestr(now, 'YYYY-MM-DD_hh-mm') '_data.mat'])
