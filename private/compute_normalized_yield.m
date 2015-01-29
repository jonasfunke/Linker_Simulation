function [ yield, p_contact] = compute_normalized_yield(ke, c_linker, d_mean, d_sigma, linker_mean, linker_sigma )
% Compute yield for varying d_mean

    

    p_contact = zeros(size(d_mean));
    
    for i=1:length(d_mean)
        p_contact(i) = colocalization_probability( d_mean(i), d_sigma, linker_mean, linker_sigma );
    end
    
    % normalize k0
    k0 = 1212/p_contact(1); % 10 bp 
    
    ki = k0 .* p_contact;
    yield = ki ./ (ki + ke.*c_linker);
    
end
