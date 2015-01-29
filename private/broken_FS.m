function [ dy ] = broken_FS( t, y, k )
% Reaction of FS with bis-maleimide linker, considering some FS are broken
% k(1) = binding from solution
%   1 = Linker, 2 = Dead, 3 = Alive, 4 = Dead+L, 5 = Alive+L  
    dy = zeros(5, 1);
    dy(1) = - k(1) * y(1) * ( y(2) + y(3) ); % Linker can bind to Dead or alive structures
    dy(2) = - k(1) * y(1) * y(2); % Dead reacts to Dead+L
    dy(3) = - k(1) * y(1) * y(3); % Alive reacts to Alive +L
    dy(4) = + k(1) * y(1) * y(2); % Dead+L
    dy(5) = + k(1) * y(1) * y(3); % Alive+L
end