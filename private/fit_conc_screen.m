%% startup
close all, clear all, clc

load('/Users/jonasfunke/Documents/Typhoon_images/2014-12-13_concentration-screen_BMH_45h/combined/2014-12-13_concentration-screen_BMH_45h_data.mat')
cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private')
c1 = [0	25	50	75	100	150	200	250	500	750	1000	5000	10000	50000	100000]/1e3;
c2 = [500	750	1000	5000	10000	50000	100000	500000	1000000	2000000	4000000	6000000	8000000	10000000]/1e3;

%%
close all
semilogx(conc_yield(1:15,1), conc_yield(1:15,2), '.', conc_yield(16:end,1), conc_yield(16:end,2), '.')

%%
i_fit = size(conc_yield,1)-6;
ft = fittype( @(k1, k2, a, x) a.*simulate_yield_vary_conc(x, [k1 k2], 0.05, 45));
fit_conc = fit(conc_yield(1:i_fit,1), conc_yield(1:i_fit,2), ft, 'startpoint',  [0.1, 100, 0.3], 'Lower', [0 0 0]) 

ci = confint(fit_conc);

dk1 = diff(ci(:,1))/2;
dk2 = diff(ci(:,2))/2;
da = diff(ci(:,3))/2;


%%
cplot = [0:0.01:0.1 0.2:0.1:1 2:10:100 200:500:5000]';
yplot = fit_conc(cplot)*100;
%%
close all
cur_fig = figure();


semilogx(cplot, yplot, '-',  conc_yield(:,1), 100.*conc_yield(:,2), 'r.', conc_yield(1:i_fit,1), 100.*conc_yield(1:i_fit,2), 'k.', ...
    cplot, 100.*fit_conc.a.*fit_conc.k2./(fit_conc.k2+cplot.*fit_conc.k1), '--'), hold on


%plot(cplot, yplot, '-',  conc_yield(:,1), 100.*(conc_yield(:,2)-conc_yield(1,2)), 'r.', conc_yield(1:i_fit,1), 100.*(conc_yield(1:i_fit,2)-conc_yield(1,2)), 'k.')
%set(gca, 'XLim', [0 1])

%cplot, 100*fit_conc.a.*fit_conc.k2./(fit_conc.k2+fit_conc.k1.*cplot), '--', ...
title({['ke = (' num2str(fit_conc.k1) '+-' num2str(dk1) ') /uM/h'],  ['ki = (' num2str(fit_conc.k2) '+-' num2str(dk2) ') /h'], ['p_intact = (' num2str(fit_conc.a*100) '+-' num2str(da*100) ') %']})

legend({'Fit to kinetic model', 'Analytical approximation', 'Data'}, 'Location', 'SE')
xlabel('Concentration [uM]')
ylabel('Yield of reaction [%]')

print(cur_fig, '-dtiff', '-r500', [path_out filesep 'concentration_screen_simple_model_01.tif'])


%%


