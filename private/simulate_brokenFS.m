function [ yield ] = simulate_brokenFS( k1, clinker, cFS, f, t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    yield = zeros(size(t));

    for i=1:length(t)
        if t(i) > 0
            time_limits = [0, t(i)];

            c_init = [clinker(i), cFS*f, cFS*(1-f), 0, 0];
            
            
        
            tmp_function = @ (t, y) broken_FS(t, y, k1);
            [T, C] = ode15s(tmp_function, time_limits, c_init);

            yield(i) = C(end, 5) ./ ( C(end,2) + C(end,3) + C(end,4) + C(end,5));

        %    close all
        %    plot(T, C(:,1), '.-', T, C(:,2), '.-', T, C(:,3), '.-', T, C(:,4), '.-', T, C(:,5), '.-')
        %    legend({'Linker', 'Dead', 'Alive', 'DL', 'AL'})
        %    set(gca, 'YLim', [0 1.1*max(c_init(2:3))], 'XLim', time_limits)
        else
            yield(i) = 0;
        end
    end
    
end

