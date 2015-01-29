%% Plots for Group-meeting 2014-12-08
%% startup
close all, clear all, clc
path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/figures';
mkdir(path_out)

%%
ke = 1.8; % /uM/h
ki = 1212; % /h
c_init = [1 0.05 0 0 0 ]; %uM
t_max = 24;

[ yield, T, C ] = simulate_yield( c_init, [ke ki], t_max, 0 );


T = T;
C = C.* 1000;
close all
plot(T, C(:,1), '.-', T, C(:,2), '.-', T, C(:,3), '.-', T, C(:,4), '.-', T, C(:,5), '.-')
legend({'Linker', 'F', 'FL', 'FL_2', 'F_LL'})
set(gca, 'YLim', [0 60], 'XLim', [0 5])
xlabel('Time [h]')
ylabel('Concentration [nM]')



%% equal rate s01
ke = 10; % /uM/h
ki = 10; % /h
c_init = [0.05 0.05 0 0 0 ]; %uM
t_max = 24;

[ yield, T, C ] = simulate_yield( c_init, [ke ki], t_max, 0 );


T = T;
C = C.* 1000;
close all
cur_fig = figure();
plot(T, C(:,1), '.-', T, C(:,2), '.-', T, C(:,3), '.-', T, C(:,4), '.-', T, C(:,5), '.-')
legend({'Linker', 'F', 'FL', 'FL_2', 'F_LL'})
set(gca, 'YLim', [0 60], 'XLim', [0 5])
xlabel('Time [h]')
ylabel('Concentration [nM]')
title({['k_e = ' num2str(ke) '/uM/h , k_i = ' num2str(ki) '/h'], ['Linker_init = ' num2str(c_init(1)) 'uM' ]})

print(cur_fig, '-dtiff', '-r500', [path_out filesep 'equal_rates_01.tif'])

%% equal rates 02
ke = 10; % /uM/h
ki = 10; % /h
c_init = [1 0.05 0 0 0 ]; %uM
t_max = 24;

[ yield, T, C ] = simulate_yield( c_init, [ke ki], t_max, 0 );


T = T;
C = C.* 1000;
close all
cur_fig = figure();
plot(T, C(:,1), '.-', T, C(:,2), '.-', T, C(:,3), '.-', T, C(:,4), '.-', T, C(:,5), '.-')
legend({'Linker', 'F', 'FL', 'FL_2', 'F_LL'})
set(gca, 'YLim', [0 60], 'XLim', [0 0.5])
xlabel('Time [h]')
ylabel('Concentration [nM]')
title({['k_e = ' num2str(ke) '/uM/h , k_i = ' num2str(ki) '/h'], ['Linker_init = ' num2str(c_init(1)) 'uM' ]})

print(cur_fig, '-dtiff', '-r500', [path_out filesep 'equal_rates_02.tif'])

