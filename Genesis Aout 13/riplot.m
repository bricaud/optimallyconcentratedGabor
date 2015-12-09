function [ ] = riplot( x )
% RIPLOT:   plot real and imaginary parts
%   usage:  riplot(x)

plot(real(x));
hold on;
plot(imag(x),'r');
hold off

end

