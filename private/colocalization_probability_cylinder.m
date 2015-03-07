function [ I ] = colocalization_probability_cylinder( d_mean, d_sigma, linker_mean, linker_sigma )
%Computes the colocalization probability of thiol and maleimdide in
%cylinder coordinates

    lim = max(d_mean+3*d_sigma, linker_mean+3*linker_sigma);
    %disp(['Limit is ' num2str(lim) ' nm.'])
    
    limits = [0, sqrt(2*lim^2), 0, lim];
    
    % pdf of thiol-thiol distance vector
    p_thiol_distance = @(a, z) (2.*pi.*a/(sqrt(2.*pi)^3.*d_sigma^3)).*exp(-a.^2/d_sigma^2/2).*exp(-(z-d_mean).^2/d_sigma^2/2);
    
    % pdf of linker
    end_to_end = @(r) pdf('Normal' ,r , linker_mean, linker_sigma);
    N = integral(@(r) 4*pi*r.^2.*end_to_end(r) ,0 , Inf); % normalization 
    p_maleimide = @(a,z) end_to_end(sqrt(a.^2+z.^2))./N;
    
    % calculate overlap
    coloc_prob = @(a,z) p_thiol_distance(a,z).*p_maleimide(a,z).*1e10;
    I = integral2(coloc_prob, limits(1), limits(2), limits(3), limits(4), 'AbsTol', 1e-24).*1e-10;  % use adaptive simpson to integrate
    
    %close(cur_fig)
    
end