function [ yield, p_contact] = compute_normalized_yield_hist_03( d_mean, d_sigma, p0,linker_parameters )
% Compute yield for varying d_mean

    p_contact = zeros(size(d_mean));
        
    for i=1:length(d_mean)
        p_contact(i) = colocalization_probability_cylinder_hist( d_mean(i), d_sigma, linker_parameters);
    end
    
    yield = 1 ./ (1 + p0./p_contact);
end
