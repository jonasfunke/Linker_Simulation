function [ dy ] = cross_linking(t, y, k)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    dy = zeros(4, 1);
    dy(1) = - k(1) * y(1) * y(2) - k(2) * y(3) * y(1);   % A = thiol-oligo
    dy(2) = - k(1) * y(1) * y(2);                        % B = maleimide-x-linker
    dy(3) = + k(1) * y(1) * y(2) - k(2) * y(3) * y(1);   % AB = oligo+linker
    dy(4) = + k(2) * y(3) * y(1);   % A2B = oligo+linker+oligo
end

