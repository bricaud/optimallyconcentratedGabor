
NN=12000;
s=0.1; % shear or slope of the chirp
e=18; % excentricity of the gaussian
      % (length/width ratio)
y = chirplet (NN, s, e); %generate the chirp
z=exp(1i*2*pi*NN/8*(1:NN)/NN); %for frequency shift
%y=y.*z; %shift to the positive frequency
cc=10; %curvature of the chirp
y=y.*exp(1i*2*exp(cc*(1:NN)/NN));


for i=1:1:2000
    cc=chirplet(12000,(i-1)*0.001,d);
   % yy=tfcenter(y);
    M(i)=cc*yy.';
end
figure()
plot(abs(M))
