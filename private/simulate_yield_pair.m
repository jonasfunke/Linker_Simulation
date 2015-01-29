function [ sim_yield ] = simulate_yield_pair( rates, c0_FS, c0_linker, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    sim_yield = zeros(size(c0_linker));
       
    for i=1:length(c0_linker)
        sim_yield(i) = simulate_yield([c0_linker(i), c0_FS, 0, 0, 0] , rates, t(i), 0);
    end

end

