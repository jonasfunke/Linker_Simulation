function [ yield, p_contact] = compute_yield(k0, ke, c_linker, d_mean, d_sigma, linker_mean, linker_sigma )
% Compute yield for varying d_mean

    yield = zeros(size(d_mean));
    p_contact = zeros(size(d_mean));

    for i=1:length(d_mean)
        p_contact(i) = colocalization_probability( d_mean, d_sigma, linker_mean, linker_sigma );
        ki = k0 .* p_contact(i);
        yield(i) = ki ./ (ki + ke.*c_linker);
    end
end

