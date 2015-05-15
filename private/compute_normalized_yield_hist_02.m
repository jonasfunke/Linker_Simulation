function [ yield, p_contact] = compute_normalized_yield_hist_02(ke, c_linker, d_mean, d_sigma, linker_parameters )
% Compute yield for varying d_mean

    p_contact = zeros(size(d_mean));
        
    for i=1:length(d_mean)
        p_contact(i) = colocalization_probability_cylinder_hist( d_mean(i), d_sigma, linker_parameters);
    end
    
    % normalize k0
    k0 = 1212/p_contact(1); % 10 bp 
    
    ki = k0 .* p_contact;
    yield = ki ./ (ki + ke.*c_linker);
    
end
