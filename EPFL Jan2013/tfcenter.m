function [ y ] = tfcenter( x )
% TFCENTER: center in time and frequency
%   usage: y = tfcenter(x)
% The calculation is based on the formula for the mean value on a circle
% mean=angle(sum_t exp(2*1i*pi*t/L).*abs(x(t))^2)
% for the mean in frequency, replace x by its fourier transform

L = length(x);

tt = 0:(L-1);
tt = exp(2*1i*pi*tt/L).';

c=1;
if ~iscolumn(x) %circshift shifts only column vectors
    x = x.';
    c=0;
end

%t_offset = round(sum(tt.*abs(x).^2)/sum(abs(x).^2));
x = x/norm(x);
%t_offset = round(sum(tt.*abs(x).^2));
t_offset = round(angle(sum(tt.*abs(x).^2))*L/2/pi);

y = circshift(x,t_offset);

ychap = fft(y);
f_offset = round(angle(sum(tt.*abs(ychap).^2)/sum(abs(ychap).^2))*L/2/pi);

ychap = circshift(ychap,f_offset);
y = ifft(ychap);

if c==0
    y = y.';
end
