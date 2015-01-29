function [ kappa, I ] = reaction_rate( k_internal, point_A, point_B, limits , plot_dist)
% k_internal = internal rate of thiol+maleimide reaction
% point_A = function handle to give  prob. dist. of point A
    
%{
    coloc_prob = @(x) point_A(x).*point_B(x);
    I = integral(coloc_prob, limits(1), limits(2));  % use adaptive simpson to integrate
    kappa = k_internal .* I;
       
    if plot_dist
        figure
        xplot = limits(1):(limits(2)-limits(1))/1000:limits(2);
        plot(xplot, point_A(xplot), xplot, point_B(xplot), xplot, sqrt(coloc_prob(xplot))  )
        legend({'point_A', 'point_B', 'sqrt(Colocalization Prob.)'})
        title(['Integral is ' num2str(I)])
    end


%}
    

    
    coloc_prob = @(x,y,z) point_A(x,y,z).*point_B(x,y,z);
    I = integral3(coloc_prob, limits(1), limits(2), limits(3), limits(4), limits(5), limits(6));  % use adaptive simpson to integrate
    kappa = k_internal .* I;
       

end

