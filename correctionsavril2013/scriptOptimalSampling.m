i = 1;
LpadV = [];
aV = [];
MV = [];
LtV = [];
xR = [1 : 0.05 : 12];
for R = 1 : 0.05 : 12
    [Lpad a M Lt] = optimalsampling_new(83001, R, 100, 100);
    LpadV = [LpadV Lpad];
    aV = [aV a];
    MV = [MV M];
    LtV = [LtV   Lt(1)* (Lpad/M)/(Lt(2)*a)];
end;    

figure;
plot(xR, LpadV);
figure;
plot(xR, aV);
figure;
plot(xR, MV);
figure;
plot(xR, LtV);
