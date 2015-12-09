function [g d s] = FindClosestChirpedDilatedGaussian2(obj)
% Find the closest generalized gaussian window from the window given in
% parameter
% Closest means: the window which maximizes the scalar product with the
% objective window
% Optimization is done by gradient ascend
% obj: objective window to be estimated by a generalized gaussian
% g : closest generalized gaussian from the objective window
% d : square root of the excentricity of the gaussian 
% s : shear of the chirp
epsilonS = 1e-5;
epsilonD = 1e-3;
startS = 0;
startD = 10;
L = length(obj);
t = -(L-1)/2 : 1 : (L-1)/2;
t = t/sqrt(L);
g_obj = fftshift(obj.');
d = startD+2*epsilonD;
new_d = startD;
s = startS+2*epsilonS;
new_s = startS;
%Objective function
f_obj = sum(conj(g_obj) .*  chirplet(L, new_s, new_d) + g_obj .*  conj(chirplet(L, new_s, new_d)));

%Gradient Loop
while (abs(new_s - s) > epsilonS || abs(new_d - d) > epsilonD)
   
   alpha = 1000; 
   s = new_s;
   d = new_d;
      
   %Derivation according to s
   uprime = conj(g_obj) .* (1i *  (t.^2) * pi) .*  sqrt(sqrt(2/L)) * (1/sqrt(abs(d))) .* exp(-pi * (t/d).^2) .* (exp(pi * i * s * t.^2));
   uprime2 = -g_obj .* (1i *  (t.^2) * pi) .*  sqrt(sqrt(2/L)) * (1/sqrt(abs(d))) .* exp(-pi * (t/d).^2) .* (exp(-pi * i * s * t.^2));
   grad_s = ((L+1)/L) * (uprime + uprime2);%./sqrt(sum(1/sqrt(pi) * 1/abs(d) * exp(-(t/d).^2)));   
   %Gradient step
   new_s = s + alpha * sum(grad_s);
  
   %Dérivation according to d
   
   %Derivation de 1/sqrt(abs(d)) .* exp(-pi * (t/d).^2) = u_ext * v_ext
   u_ext =  1/sqrt(abs(d));
   v_ext = exp(-pi * (t/d).^2);
   u_extprime = ((-0.5*sign(d))/abs(d)^(3/2));
   v_extprime = 2*pi*exp(-pi * (t/d).^2) .* t.^2 / d^3 ;
   grad_ext_d = u_ext * v_extprime + u_extprime * v_ext;
   
   %Derivation of conj(g_obj) .* g
   grad_d_factor = conj(g_obj) .* sqrt(sqrt(2/L)) .* ((L+1) /L) .*exp(pi * i * s * t.^2);
    
   %Derivation de g_obj .* conj(g)
   grad_d_factor2 = g_obj .* sqrt(sqrt(2/L)) .* ((L+1) /L) .* exp(-pi * i * s * t.^2);
   
   %Gradient according d
   grad_d = (grad_d_factor + grad_d_factor2) .* grad_ext_d;

   %Update of excentricity/dilation
   new_d = d + alpha * sum(grad_d);
   
   %Objective function
   newf_obj = 	sum(conj(g_obj) .*  chirplet(L, new_s, new_d) + g_obj .*  conj(chirplet(L, new_s, new_d)));

   %Looking for a rather big alpha to progress fast enough
   %while we do not progress, then divide alpha by 2
   while (newf_obj <= f_obj && alpha > 1e-8)
      alpha = alpha/2;
      new_d = d + alpha * sum(grad_d);    
      new_s = s + alpha * sum(grad_s);
      newf_obj = sum(conj(g_obj) .*  chirplet(L, new_s, new_d) + g_obj .*  conj(chirplet(L, new_s, new_d)))   ;
   end;
   if newf_obj > f_obj
       f_obj = newf_obj;
   else
      new_d = d;
      new_s = s;
   end

   %Try to progress only on the direction of the gradient of d
   alpha = 1000;
   d3 = new_d;
   while (newf_obj <= f_obj && alpha > 1e-8)
      alpha = alpha/2;
      new_d = d3 + alpha * sign(sum(grad_d));
      newf_obj = sum(conj(g_obj) .*  chirplet(L, new_s, new_d) + g_obj .*  conj(chirplet(L, new_s, new_d)))   ;
   end;
   if newf_obj > f_obj
       f_obj = newf_obj;
   else
      new_d = d3;
   end 
   
   
   fprintf('Fmin: %f\tS: %f\tD: %f\n', f_obj, new_s, new_d)
   
   %Try to progress only on the direction of the gradient of s
   alpha = 1000;
   s2 = new_s;
   while (newf_obj <= f_obj && alpha > 1e-8)
      alpha = alpha/2;
      new_s = s2 + alpha * sign(sum(grad_s));
      newf_obj = sum(conj(g_obj) .*  chirplet(L, new_s, new_d) + g_obj .*  conj(chirplet(L, new_s, new_d)))   ;
   end;
   if newf_obj > f_obj
       f_obj = newf_obj;
   else
      new_s = s2;
   end
    

end
%Construct the window to be returned
g = fftshift(chirplet(L, s, d)).';