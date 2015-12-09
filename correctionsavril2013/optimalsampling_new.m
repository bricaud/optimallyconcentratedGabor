function [Lpad,a,M,lt] = optimalsampling_new(L,R,e,sh)
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
Ac(2,1)=Ac(2,1)-floor(Ac(2,1)/Ac(2,2))*Ac(2,2);

ac = Ac(1,1);
bc = Ac(2,2);
sc = Ac(2,1);
sc = mod(sc,bc);
ltc = sc/bc;
A = Ac*sqrt(L);

mg=1.05*L;
L=nextfastfft(L);
while L(end)<mg
	L=[L, nextfastfft(L(end)+1)];
end


nlen=length(L);
for ii=1:nlen
	ldiv=extfactor(L(ii));
	[~,inda]=min(abs(ldiv-A(1,1)));
	[~,indb]=min(abs(ldiv-A(2,2)));
	a(ii)=ldiv(inda);
	b(ii)=ldiv(indb);
	mins=a(ii)*b(ii)/gcd(a(ii)*b(ii),L(ii));
	s(ii)=ceil((A(2,1)-mins/2)/mins)*mins;
	err(ii)=norm(A-[a(ii) 0;s(ii) b(ii)]);
end

[~,opt]=min(err);

Lpad=L(opt);
a=a(opt);
b=b(opt);
s=s(opt);
lt=[s, b]/gcd(b,s);
M=Lpad/b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ldiv pfacs] = extfactor(L)

tmp=factor(L);
pfacs=unique(tmp);
ll=length(pfacs);
pfacs=[pfacs;zeros(1,ll)];

for k=1:ll
	pfacs(2,k) = sum(tmp==pfacs(1,k));
end

ldiv=pfacs(1,1).^(0:pfacs(2,1));
k=1;
for f=pfacs(1,2:end)
    k=k+1;
	ldiv=ldiv'*(f.^(0:pfacs(2,k)));
	ldiv=ldiv(:)';
end
% ldiv=sort(ldiv);
