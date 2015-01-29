%%
clear all
close all
clc
path_out = '/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/figures';

cd('/Users/jonasfunke/Dropbox/POSITIONING/linker_simulation/private')
%% set parameters
%d = 3; % nm, distance set by spacer length
sigma_thiol = 0.5; % nm, fluctuations of thiol and arm of FS

L = [8.18 10.16 10.65 11.13 14.53]/10 ;   % nm, BMOE, BMH mPDM, pPDM, BM
sigma_L = [0.75 2.41 0.55 0.52 1.51] / 10; % nm, fluctuations of maleimide linker
linker_names = {'BMOE', 'BMH', 'mPDM', 'pPDM', 'BM'};

p_react = 1; % reaction rate if hit occurs

%
n = [0:50]';
d = n.*0.335*16/47;

rate2 = zeros(length(d), length(L));

close all
figure
disp('Simulation started')
for j=1:length(L)
    disp(['iteration ' num2str(j) ' of ' num2str(length(L))])
for i=1:length(d)
    %disp(['iteration ' num2str(i) ' of ' num2str(length(d))])

    % prob dist of thiol
    thiol = @(x, y, z) pdf('Normal',x,0,sigma_thiol).*pdf('Normal',y,0,sigma_thiol).*pdf('Normal',z,d(i),sigma_thiol);

    %prob dist of maleimide linker
    end_to_end = @(r) pdf('Normal',r,L(j),sigma_L(j));
    N = integral(@(r) 4*pi*r.^2.*end_to_end(r) ,0, Inf); % normalization (remember raumwinkel integration)
    maleimide = @(x,y,z) end_to_end(sqrt(x.^2+y.^2+z.^2))./N;
    
    lim(1) = L(j) - 5;
    lim(2) = L(j) + 5;
    %calcluate rate
    rate2(i,j) = reaction_rate(p_react, thiol, maleimide, [lim(1) lim(2) lim(1) lim(2) lim(1) lim(2)], 1);

    plot(d, rate2, '.-')
    pause(1)
end
end
disp('Simulation done')

%% plot coloc prob
myleg =  cell(length(L),1);

for i=1:length(L)
    myleg{i} = [linker_names{i} ' l_0 = (' num2str(L(i)) '\pm' num2str(sigma_L(i)) ') nm'];
end

close all
cur_fig = figure();
plot(d, rate2, '-'), hold on
legend(myleg )
xlabel('Separation between thiols [nm]')
ylabel('Probability of colocalization')
print(cur_fig, '-dtiff', '-r300', [path_out filesep 'internal_rate_sigma-linker-screen.tif'])

%% normalize to BMH linker at 10bp spacer
p_hit = 1212./rate2(11,2);

ki = zeros(size(rate2));
for i=1:size(rate2,2)
    ki(:,i) = rate2(:,i).*p_hit;
end

close all
cur_fig = figure();
subplot(2,1,1)
plot(d, ki, '-'), hold on
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
plot(d, yield, '-'), hold on
xlabel('Separation between thiols [nm]')
ylabel('Simulated yield')
legend(myleg )
set(gca, 'YLim', [0 1.2])

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'internal_rate_sigma-linker-screen02.tif'])

%%

close all
ke = 0.18;
c_linker = [1 10 50 100 1000]; %uM

yield = zeros(size(rate2));
for i=1:size(yield,2)
    subplot(5,1,i)
    for j=1:5
        yield = ki(:,i)./ (ki(:,i) + c_linker(j).*ke);
        plot(d, yield, '-'), hold on
    end
end
xlabel('Separation between thiols [nm]')
ylabel('Simulated yield')
set(gca, 'YLim', [0 1.2])








