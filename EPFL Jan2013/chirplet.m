%Create a dilated chirped gaussian, centered in L/2
function g = chirplet (L, s, e)
% L: Length of the window
% s: Frequency shear of the chirp
% e: time dilation/excentricity
% Time vector 
t = -(L-1)/2 : 1 : (L-1)/2;
t = t/sqrt(L);
%Create a dilated gaussian
g = sqrt(sqrt(2/L)) * 1/sqrt(abs(e)) * (exp(-pi * t.^2 / e^2)) ;
% Chirp the gaussian
g = (L+1)/L * g .* (exp(1i * pi * s * (t.^2)));
