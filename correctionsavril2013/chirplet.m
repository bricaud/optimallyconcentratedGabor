%Create a dilated chirped gaussian, centered in L/2
function g = chirplet (L, s, e)
% L: Length of the window
% s: Frequency shear of the chirp
% e: time dilation/excentricity
% Time vector 
t = -(L)/2 : 1 : (L-1)/2;
t = t/sqrt(L);
% no negative variance accepted:
e=abs(e);
%Create a dilated gaussian, compatible with pgauss of LTFAT
g = sqrt(sqrt(2/L)) * 1/sqrt(sqrt(e)) * (exp(-pi * t.^2 / e)) ;
% Chirp the gaussian
g = g .* (exp(1i * pi * s * (t.^2)));
