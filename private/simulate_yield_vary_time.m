function [ sim_yield ] = simulate_yield_vary_time(rates, c_init, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    sim_yield = zeros(length(t),1);
    for i=1:length(t)
        sim_yield(i) = simulate_yield(c_init , rates, t(i), 0);
    end

end

