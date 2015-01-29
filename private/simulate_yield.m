function [ yield, T, C ] = simulate_yield( c_init, rates, t_max, plot_simulation )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%c_inint =  Linker, FS, FSL_open, FSL2, FSL_closed in uM
%rates =  on-Rate from solution [1/h/uM], internal closing rate [1/h]


if t_max > 0
    time_limits = [0 t_max]; % start and end time in seconds


    %options = odeset('RelTol',1e-4);%,'AbsTol',1e-1 .* [1 1 1 1 1]);
    %options = odeset('AbsTol',[1 1 1 1 1]);

    tmp_function = @ (t, y) linking_reaction(t, y, rates);

    [T, C] = ode15s(tmp_function, time_limits, c_init);

    yield = C(end, 5) ./ ( C(end,2) + C(end,3) + C(end,4) + C(end,5));

    time_99 = T(diff(C(:,4)+C(:,5) > 0.99*min(c_init(1:2)))==1); % get time when 99% o reaction is finished
    if isempty(time_99)
        time_plot = time_limits;
     %   disp('Warning: Simulation not converged!')
    else
        time_plot = [0 time_99];
      %  disp(['99 % reached after ' num2str(time_plot(2)) ' s'])

    end
    if plot_simulation
        figure;
        plot(T, C(:,1), '.-', T, C(:,2), '.-', T, C(:,3), '.-', T, C(:,4), '.-', T, C(:,5), '.-')
        legend({'Linker', 'F', 'FL', 'FL2', 'FL_linked'})
        set(gca, 'YLim', [0 1.1*c_init(2)], 'XLim', time_plot)
        xlabel('Time [h]')
        ylabel('Concentration [uM]')
    end
else
    
    yield = 0;
    T = zeros(5,0);
    C = zeros(5,0);
end

end

