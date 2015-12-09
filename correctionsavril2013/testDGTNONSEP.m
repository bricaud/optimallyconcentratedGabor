[yfull FS] = wavread('RPMRaiseSquare10Harm.wav');
y = yfull(20001:160000);
L = length(y);
Lpad = 147000;
f = [y ; zeros(Lpad - L, 1)];
%Desired redundancy of the dgt transform
R = 32;
%initial window
d = 716;
s = 0.0045;
a = 1225; % time step
M = 36750; % number of frequency channels
gamma = chirplet(Lpad, s, d); % Example of window

for i = 1: 1
    i, 1
    tmp1 = abs(dgt(gamma, f, a, M));
    i, 2
        tmp1 = abs(dgt(gamma, f, a, M, 'lt', [4,7]));
        i, 3
        tmp1 = abs(dgt(gamma, f, a, M, 'lt', [1,4]));
        i, 4
        tmp1 = abs(dgt(gamma, f, a, M, 'lt', [3,4]));
        i, 5
        tmp1 = abs(dgt(gamma, f, a, M, 'lt', [1,10]));
        i, 6
        tmp1 = abs(dgt(gamma, f, a, M, 'lt', [0,1])); 
end;