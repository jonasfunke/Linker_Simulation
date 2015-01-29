function [ dy ] = linking_reaction( t, y, k )
% Reaction of FS with bis-maleimide linker
% k(1) = binding from solution
% k(2) = binding if on end is already attached
%   1 = Linker, 2 = FS, 3 = FS with 1 Linker - open, 4 = FS with 2 linker - open, 5 = FS with 1 linker-closed  
    dy = zeros(5, 1);
    dy(1) = - k(1) * y(1) * ( y(2) + y(3)); % Linker can bind to FS and to FSL-open
    dy(2) = - k(1) * y(1) * y(2); % FS only
    dy(3) = + k(1) * y(1) * y(2) - k(1) * y(3) * y(1) - k(2) * y(3); % FS with 1 linker, open
    dy(4) = + k(1) * y(3) * y(1); % FSL2 - open
    dy(5) = + k(2) * y(3); % FSL, closed
end

