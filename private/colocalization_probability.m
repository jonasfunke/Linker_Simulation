function [ I ] = colocalization_probability( d_mean, d_sigma, linker_mean, linker_sigma )
%Computes the colocalization probability of thiol and maleimdide

    lim = max(d_mean+3*d_sigma, linker_mean+3*linker_sigma);
    %disp(['Limit is ' num2str(lim) ' nm.'])
    
    limits = [-lim, lim, -lim, lim, 0, lim];
    
    % pdf of thiol-thiol distance vector
    p_thiol_distance = @(x, y, z) pdf('Normal',x,0,d_sigma).*pdf('Normal',y,0,d_sigma).*pdf('Normal',z,d_mean,d_sigma);
    
    % pdf of linker
    end_to_end = @(r) pdf('Normal' ,r , linker_mean, linker_sigma);
    N = integral(@(r) 4*pi*r.^2.*end_to_end(r) ,0 , Inf); % normalization 
    p_maleimide = @(x,y,z) end_to_end(sqrt(x.^2+y.^2+z.^2))./N;

    % plot along z-axis
    %{
    cur_fig = figure();
    zplot = 0:lim/1000:lim;
    plot(zplot, p_thiol_distance(0, 0, zplot), zplot, end_to_end(zplot)./N )
    legend({'Free thiol', 'Free maleimide'})
    xlabel('Distance between thiols (z-axis direction)')
    ylabel('Probability')
    %}
    
    % calculate overlap
    coloc_prob = @(x,y,z) p_thiol_distance(x,y,z).*p_maleimide(x,y,z);
    I = integral3(coloc_prob, limits(1), limits(2), limits(3), limits(4), limits(5), limits(6));  % use adaptive simpson to integrate
    
    %close(cur_fig)
    
end

