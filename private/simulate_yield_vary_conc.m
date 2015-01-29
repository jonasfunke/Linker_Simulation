function [ sim_yield ] = simulate_yield_vary_conc( c_linker, rates, c_FS, t_max )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    c_init = [c_linker c_FS*ones(size(c_linker,1),1) zeros(size(c_linker,1),1) zeros(size(c_linker,1),1) zeros(size(c_linker,1),1)];
    sim_yield = zeros(size(c_linker,1),1);
    for i=1:size(c_init,1)
        sim_yield(i) = simulate_yield(c_init(i,:) , rates, t_max,0);
    end

end

