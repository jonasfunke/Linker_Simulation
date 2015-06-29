function [ prob_density ] = linker_end_to_end( r, p )
%uses the end-to-end distance distribution to compute the real  3d
%locaizatzion probability density
%   Detailed explanation goes here

    N = integral(@(x) linker_end_to_end_dist(x, p) ,0 , Inf); % normalization int 4 pi r2 dr  p_endtoend / (4 pi r2) 
    prob_density = linker_end_to_end_dist( r, p )./N./(4.*pi.*r.^2);
    
    
end

