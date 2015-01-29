%%
clear all
close all
clc
path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/figures';

cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private')
%% set parameters
%d = 3; % nm, distance set by spacer length
sigma_thiol = [0.01 0.1 0.2 0.3 0.4 0.5 0.75 1 2]; % nm, fluctuations of thiol and arm of FS

L = 1.016;   % nm, length of BMH linker
sigma_L = 0.24; % nm, fluctuations of maleimide linker

p_react = 1; % reaction rate if hit occurs

%
n = [0:50]';
d = n.*0.335*16/47;

rate2 = zeros(length(d), length(sigma_thiol));

close all
figure
disp('Simulation started')
for j=1:length(sigma_thiol)
    disp(['iteration ' num2str(j) ' of ' num2str(length(sigma_thiol))])
for i=1:length(d)
    %disp(['iteration ' num2str(i) ' of ' num2str(length(d))])

    % prob dist of thiol
    thiol = @(x, y, z) pdf('Normal',x,0,sigma_thiol(j)).*pdf('Normal',y,0,sigma_thiol(j)).*pdf('Normal',z,d(i),sigma_thiol(j));

    %prob dist of maleimide linker
    end_to_end = @(r) pdf('Normal',r,L,sigma_L);
    N = integral(@(r) 4*pi*r.^2.*end_to_end(r) ,0, Inf); % normalization (remember raumwinkel integration)
    maleimide = @(x,y,z) end_to_end(sqrt(x.^2+y.^2+z.^2))./N;
    
    lim(1) = L - 5;
    lim(2) = L + 5;
    %calcluate rate
    rate2(i,j) = reaction_rate(p_react, thiol, maleimide, [lim(1) lim(2) lim(1) lim(2) lim(1) lim(2)], 1);

    plot(d, rate2, '.-')
    pause(1)
end
end
disp('Simulation done')

%% plot coloc prob
myleg =  cellstr(num2str(sigma_thiol(:)));
for i=1:length(sigma_thiol)
    myleg{i} = ['\sigma_{thiol} = ' myleg{i} ' nm'];
end

close all
cur_fig = figure();
plot(d, rate2, '-'), hold on
legend(myleg )
xlabel('Separation between thiols [nm]')
ylabel('Probability of colocalization')
print(cur_fig, '-dtiff', '-r300', [path_out filesep 'internal_rate_sigma-thiol-screen.tif'])

%% plot coloc prob
myleg =  cellstr(num2str(sigma_thiol(:)));
for i=1:length(sigma_thiol)
    myleg{i} = ['\sigma_{thiol} = ' myleg{i} ' nm'];
end
myleg = myleg(2:end);
p_hit = 1212./rate2(11,:);

ki = zeros(size(rate2));
for i=1:size(rate2,2)
    ki(:,i) = rate2(:,i).*p_hit(i);
end

close all
cur_fig = figure();
subplot(2,1,1)
plot(d, ki(:,2:end), '-'), hold on
vline(d(11), 'k--');
legend(myleg )
xlabel('Separation between thiols [nm]')
ylabel('Estimated internal Rate [1/h]')

subplot(2,1,2)
ke = 0.18;
c_linker = 10; %uM

yield = zeros(size(rate2));
for i=1:size(yield,2)
    yield(:,i) = ki(:,i)./ (ki(:,i) + c_linker.*ke);
end
plot(d, yield(:,2:end), '-'), hold on
xlabel('Separation between thiols [nm]')
ylabel('Internal Rate [1/h]')
legend(myleg )
set(gca, 'YLim', [0 1.2])

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'internal_rate_sigma-thiol-screen02.tif'])


%%

ki = [0:1500]'; %/h
ke = 0.18; % /uM/h
c = [1 10 50 100 1000 2000]; %uM
yield = zeros(length(ki), length(c));

for i=1:length(c)
    yield(:,i) = ki./(ki+c(i).*ke);
end
myleg =  cellstr(num2str(c'));
for i=1:length(c)
    myleg{i} = ['c_{linker} = ' myleg{i} ' uM'];
end

close all
cur_fig = figure();
plot(ki, yield)
legend(myleg, 'location', 'SE')

xlabel('Internal Rate k_i [1/h]')
ylabel('Yield of reaction')
title(['k_e = ' num2str(ke) '/uM/h' ])
print(cur_fig, '-dtiff', '-r300', [path_out filesep 'kinetic_model_approximation01.tif'])

%%
close all
yield = zeros(length(rate2),1);
for i=1:length(rate2)
    yield(i) = simulate_yield(c_init, [1e-4 rate2(i)], 0);
end

plot(d,yield)
%%


conc_screen =[    0.3700    0.3054
    0.0412    0.0137
    0.0787    0.0491
    0.1386    0.1044
    0.1721    0.1376
    0.2074    0.1697
    0.2860    0.2303
    0.3295    0.2743
    0.3489    0.2945
    0.4159    0.3530
    0.4272    0.3564
    0.4275    0.3606
    0.3997    0.3425
    0.3197    0.2773
    0.4033    0.3607];

 data = [1e-3.*[0 25 50 75 100 150 200 250 500 750 1e3 10e3 2e6]' conc_screen(2:end-1,2)];

 close all
cur_fig = figure();

semilogx(data(:,1), data(:,2).*100, '.-')

xlabel('Concentration [uM]')
ylabel('Yield of reaction [%]')

print(cur_fig, '-dtiff', '-r500', [path_out filesep 'concentration_screen_data-only_01.tif'])

 
 

%%
ft = fittype( @(k1, k2, a, x) a.*simulate_yield_vary_conc(x, [k1 k2], 0.05, 46) );
fit_conc = fit(data(:,1), data(:,2), ft, 'startpoint',  [0.1, 1e2, 0.4], 'Lower', [0 0 0]) 

ci = confint(fit_conc);

dk1 = diff(ci(:,1))/2;
dk2 = diff(ci(:,2))/2;
da = diff(ci(:,3))/2;

cplot = [0:0.01:0.1 0.2:0.1:1 2:10:100 200:500:5000]';

%%
close all
cur_fig = figure();

semilogx(cplot, fit_conc(cplot)*100, '-', ...
    cplot, 100*fit_conc.a.*fit_conc.k2./(fit_conc.k2+fit_conc.k1.*cplot), '--', ...
    data(:,1), 100*data(:,2), '.k')
title({['ke = (' num2str(fit_conc.k1) '+-' num2str(dk1) ') /uM/h'],  ['ki = (' num2str(fit_conc.k2) '+-' num2str(dk2) ') /h'], ['p_intact = (' num2str(fit_conc.a*100) '+-' num2str(da*100) ') %']})

legend({'Fit to kinetic model', 'Analytical approximation', 'Data'}, 'Location', 'SE')
xlabel('Concentration [uM]')
ylabel('Yield of reaction [%]')

print(cur_fig, '-dtiff', '-r500', [path_out filesep 'concentration_screen_simple_model_01.tif'])

%%
k1 = 0.18;
k2 = 1200;
%a = 0.356;
cplot = [0:0.01:0.1 0.2:0.1:1 2:10:100 200:500:5000]';

myfun = @(x) 1.*simulate_yield_vary_conc(x, [k1 k2], 0.05, 46) ;

close all
cur_fig = figure();
semilogx(cplot, myfun(cplot), '-', ...
    cplot, k2./(k2+k1.*cplot), '--',... 
    data(:,1), data(:,2)./0.356, '.k')

xlabel('Concentration [uM]')
ylabel('Yield of reaction')


%%
yield_100nM_BMH =[
    0       0
    2.5    0.0313
    7.833    0.0609
    27.5    0.1294
    48.333    0.1617
    73.5    0.1865
    142.5    0.2173
    14*24    0.2464]; % reaction time [h], yield cy5

%
yield_500nM_BMH =[
    0         0
    2.5    0.0794
    7.833    0.1881
    27.5    0.2746
    48.333    0.2865
    73.5    0.3058
    142.5    0.3087
    14*24    0.3097];% reaction time [h], yield cy5

t = yield_100nM_BMH(:,1);
y1 = yield_100nM_BMH(:,2);
y5 = yield_500nM_BMH(:,2);
ct_data = [0.1.*ones(size(t)), t, y1; 0.5.*ones(size(t)), t, y5];

%%
k2 = 1212;

%ft = fittype( @(k1, k2, c0_L, x) simulate_yield_vary_time([k1 k2], [c0_L 0.050 0 0 0], x) );
ft1 = fittype( @(k1,  x) simulate_yield_vary_time([k1 k2], [0.1 0.050 0 0 0], x) );
ft5 = fittype( @(k1, x) simulate_yield_vary_time([k1 k2], [0.5 0.050 0 0 0], x) );


t = yield_100nM_BMH(:,1);
y1 = yield_100nM_BMH(:,2)./0.31;
y5 = yield_500nM_BMH(:,2)./0.31;

time_series100 = fit(t, y1, ft1, 'startpoint', [0.2], 'Lower', [0 ]) 
time_series500 = fit(t, y5, ft5, 'startpoint', [0.2 ], 'Lower', [0 ]) 


tplot = t; %0:1:t(end);

close all
cur_fig = figure();

plot(  tplot, time_series100(tplot), tplot, time_series500(tplot), t, y1, '.', t, y5, '.'  )

legend({'Fit to 100 nM', 'Fit to 500 nM', '100 nM Linker', '500 nM Linker'}, 'Location', 'SE')
xlabel('Time [h]')
ylabel('Normalized yield')
set(gca, 'YLim', [0 1.1])
print(cur_fig, '-dtiff', '-r500', [path_out filesep 'time_series.tif'])

%%
%ft = fittype( @(k1, k2, c0_L, x) simulate_yield_vary_time([k1 k2], [c0_L 0.050 0 0 0], x) );
ft1 = fittype( @(k1, c0,  x) simulate_yield_vary_time([k1 k2], [c0 0.050 0 0 0], x) );
ft5 = fittype( @(k1, c0, x) simulate_yield_vary_time([k1 k2], [c0 0.050 0 0 0], x) );


t = yield_100nM_BMH(:,1);
y1 = yield_100nM_BMH(:,2)./0.31;
y5 = yield_500nM_BMH(:,2)./0.31;

time_series100 = fit(t, y1, ft1, 'startpoint', [0.2 0.1], 'Lower', [0 ]) 
time_series500 = fit(t, y5, ft5, 'startpoint', [0.2 0.5], 'Lower', [0 ]) 


tplot = t; %0:1:t(end);

close all
cur_fig = figure();

plot(  tplot, time_series100(tplot), tplot, time_series500(tplot), t, y1, '.', t, y5, '.'  )

legend({'Fit to 100 nM', 'Fit to 500 nM', '100 nM Linker', '500 nM Linker'}, 'Location', 'SE')
xlabel('Time [h]')
ylabel('Normalized yield')
set(gca, 'YLim', [0 1.1])
print(cur_fig, '-dtiff', '-r500', [path_out filesep 'time_series_fit-linker-conc.tif'])

%% simulate only  500nM with simple model


%% fit 500nM data with simple decay into two 
FS_tot = 0.05;



ft = fittype( @(k1, f, clinker, t) simulate_brokenFS(k1, clinker, FS_tot, f, t),...
    'independent', {'clinker', 't'});

i= 8;

fit100 = fit(ct_data(1:i,1:2), ct_data(1:i,3), ft, 'startpoint', [0.1 0.5], 'Lower', [0 0]) 
fit500 = fit(ct_data(9:end,1:2), ct_data(9:end,3), ft, 'startpoint', [0.1 0.5], 'Lower', [0 0]) 


close all
cur_fig = figure();

plot( ct_data(1:8,2), fit100(ct_data(1:8,1:2)), '-', ct_data(9:end,2), fit500(ct_data(9:end,1:2)), ...
    '-', ct_data(1:end,2), ct_data(1:end,3), '.k')
legend({'100 nM Linker', '500 nM Linker', 'Data'}, 'Location', 'best')
xlabel('Time [h]')
ylabel('Yield')
print(cur_fig, '-dtiff', '-r500', [path_out filesep 'time_series_fit-dead-FS.tif'])



%%
k2 = 1.2e3;

ft = fittype( @(k1, c0, x) simulate_yield_pair([k1, k2], 0.05, c0, x) ,...
    'independent', {'c0', 'x'});


kinetic_fit = fit(ct_data(:,1:2), ct_data(:,3), ft, 'startpoint', [1e-1], 'Lower', [0 ]) 


close all
plot(ct_data(1:8,2), kinetic_fit(ct_data(1:8,1:2)), '.-', ct_data(9:end,2), kinetic_fit(ct_data(9:end,1:2)), '.-', ...
    ct_data(:,2), ct_data(:,3), '.k')
legend({'100 nM Linker', '500 nM Linker', 'Data'}, 'Location', 'best')

%% fit 100nM data with simple decay into two 
FS_tot = 0.05;

ft = fittype( @(k1, f, clinker, t) simulate_brokenFS(k1, clinker, FS_tot, f, t),...
    'independent', {'clinker', 't'});

broken_fit = fit(ct_data(:,1:2), ct_data(:,3), ft, 'startpoint', [0.1 0.5], 'Lower', [0 0]) 

close all
plot(ct_data(1:8,2), broken_fit(ct_data(1:8,1:2)), '.-', ct_data(9:end,2), broken_fit(ct_data(9:end,1:2)), '.-', ...
    ct_data(:,2), ct_data(:,3), '.k')
legend({'100 nM Linker', '500 nM Linker', 'Data'}, 'Location', 'best')

%%
k1 = 0.01;
FS_tot = 0.05;
f = 0.68;
clinker1 = 0.2 .* ones(size(t));
clinker5 = 1 .* ones(size(t));

y1 = simulate_brokenFS(k1, clinker1, FS_tot, f, t);
y5 = simulate_brokenFS(k1, clinker5, FS_tot, f, t);

close all
plot(t, y1 , '.-', t, y5, '.-', ...
    ct_data(:,2), ct_data(:,3), '.k')


%%

