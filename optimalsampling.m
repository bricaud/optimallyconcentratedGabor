function [Lpad,a,M,lt] = optimalsampling(L,R,e,sh)
% function [Lpad,a,M,lt] = optimalsampling(L,R,e,sh)
% suggests a good (discrete) sampling pattern for given input parameters
%
% L  - length of the signal to be investigated
% R  - desired redundancy
% e  - excentricity of the generalized Gaussian
% sh - shear of the generalized Gaussian
%
%
% Output:
% Lpad - suggested signal length one has to pad to
% a    - discrete step parameter
% M    - number of channels
% lt   - sampling strategy (lattice type)


Ac = sqrt(1/R)*[1/sqrt(2)*3^.25 0; 1/(sqrt(2)*3^.25) sqrt(2)/3^.25];
Ac = diag(sqrt([e,1/e]))*Ac;
Ac = [1 0;sh 1]*Ac;

ac = Ac(1,1);
bc = Ac(2,2);
sc = Ac(2,1);
sc = mod(sc,bc);
ltc = sc/bc;

ltc = round(ltc*10)/10;
lt1 = ltc*10; lt2 = 10;
gc = gcd(lt1,lt2);

a = round(ac*sqrt(L));
b = round(bc*sqrt(L)/2)*2;
lt = [lt1/gc lt2/gc];
lt1 = lt(1); lt2 = lt(2);
s = b*lt1/lt2;

currerr = 100000;
for llt1=10:-1:1
    for llt2=10:-1:1
        feas = (b/llt2==round(b/llt2));
        feas = (feas && gcd(llt1,llt2)==1 && llt1<=llt2);
        if feas
            ss = b*llt1/llt2;
            err = abs(s-ss);
            if err<currerr
                currerr = err;
                lt1ch = llt1;
                lt2ch = llt2;
            end
        end
    end
end

lt = [lt1ch lt2ch];


Lmin = lt(2)*lcm(a,b);

Lpad = ceil(L/Lmin)*Lmin;

M = Lpad/b;
