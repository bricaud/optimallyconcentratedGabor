function [alpha_ton alpha_trans bruit] = fista(f, p1, q1, p2, q2, noisepond, epsilon, parametre, splitbysec,stdvalue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [alpha_ton alpha_trans bruit] = fista(f, p1, q1, p2, q2, epsilon,
%   parametre, splitbysec, stdvalue)
% 
% FUNCTION:
%   Perform the sparse regression in order to decompose of a wav signal into tonal, transient and noisy
%   components. Objective function: 0.5*|sig - to -tr|^2 + noisepond * (
%   lambda * |to|_1 + (1-lambda) |to_2|).
% 
% 
% INPUT:
%       f: input signal to separate
%       p1,q1: norm used to regularize tonal part
%       p2,q2: norm used to regularize transient part
%       noisepond: ponderation of noisy part
%       epsilon: stopping criterion of the optimization
%       parametre: Size of windows for STFT, overlap, sample frequency
%       splitbysec: Nb of splits (different ponderations) by second
%       stdvalue: standard divisor for tonal/transient estimator ponderation
%
% OUTPUT:
%       tonale   : tonal component of the signal
%       transitoire: transient component of the signal
%       bruit: noisy component of the signal
%
%
% NOTA: An appropriate weigth for the 2 regularization is estimated on the
% original signal
%  
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initializing and getting the parameters
fprintf('****************Initialisation*****************\n');

%Initialisation of the fista parameter
t=1;

taille_fenetre_ton = parametre.taille_fenetre_ton;
taille_fenetre_trans = parametre.taille_fenetre_trans;
noverlap_ratio = parametre.noverlap_ratio;
fe = parametre.fe;


%% Building of the Hanning windows
fprintf('****************Hanning windows*****************\n');
%First for the tonal part
fen_ton = hanning(taille_fenetre_ton);
%Then for the transient part
fen_trans = hanning(taille_fenetre_trans);

%Compute the sample overlap for each component 
noverlap_ton=ceil(noverlap_ratio*taille_fenetre_ton);
noverlap_trans=ceil(noverlap_ratio*taille_fenetre_trans);

%Build the strict window:
fen_ton=compute_tight_window(fen_ton,taille_fenetre_ton-noverlap_ton);
fen_trans=compute_tight_window(fen_trans,taille_fenetre_trans-noverlap_trans);


%% 0 Padding of the original signal
fprintf('*********************Build Signal**********************\n');

%Compute time steps
pas_temps_ton=taille_fenetre_ton-noverlap_ton;
pas_temps_trans=taille_fenetre_trans-noverlap_trans;

%0 Padding in order to get the right number of window
long_sig=length(f);
while(rem(length(f)-noverlap_ton,pas_temps_ton)~=0 || rem(length(f)-noverlap_trans,pas_temps_trans)~=0)
f=[f;0];
end


%% Compute the first TF diagrams
fprintf('****************Premier diagramme TF*****************\n');

%Compute the first TF diagram for both the tonal and the transient views
signal_spec_ton_init=spectrogram(f,fen_ton,noverlap_ton,taille_fenetre_ton,fe);
signal_spec_trans_init=spectrogram(f,fen_trans,noverlap_trans,taille_fenetre_trans,fe);

%Initialisation of the tonal component to the full signal and of the
%transient component to zero (this seems to make the optimisation faster)
signal_spec_ton=signal_spec_ton_init;
signal_spec_trans=zeros(size(signal_spec_trans_init));

%Current temporal signals
alpha_ton=ispecgram(signal_spec_ton,fe,fen_ton,noverlap_ton);
alpha_trans=ispecgram(signal_spec_trans,fe,fen_trans,noverlap_trans);

%Compute the stepsize of the gradient descent step of FISTA
stepsize = compute_stepsize(f,signal_spec_ton,signal_spec_trans);



%% Compute the 
%% regularization weights

%Around one split (one different regularization weight) per one half-second
nbsplits=ceil((splitbysec)*length(f)/fe);

%Get the ponderation vectors
[pond1vec pond2vec] = compute_ponderation(signal_spec_ton_init,signal_spec_trans_init,nbsplits,taille_fenetre_ton,taille_fenetre_trans,noverlap_ratio,fe,noisepond,stdvalue);



%% Initialization of optimization

%Current Spectrograms
spec_ton=signal_spec_ton;
spec_trans=signal_spec_trans;
%Spectrograms for the next step
spec_tonnext=signal_spec_ton;
spec_transnext=signal_spec_trans;
%Spectrograms for the last step
spec_tonlast=signal_spec_ton;
spec_translast=signal_spec_trans;






%% Creation of groups for the norm

%For tonal, no group, norm 1
fprintf('****************Initialisation Groupe 1*****************\n');
group1={1:numel(spec_ton)};
%For transient, groups of length 1 (1 frame), mixed norm 2,1
fprintf('****************Initialisation Groupe 2*****************\n');
group2=create_groups_bands(1,size(spec_trans,1),size(spec_trans,2),fe);


%% Ready for optimization
fprintf('****************Iterations*****************\n');

%First iteration

%Compute the objective function terms
[obj1 obj2 obj3] = compute_obj_function_signal(f,alpha_ton,alpha_trans,spec_ton,spec_trans,p1,q1,p2,q2,group1,group2);

%Number of iterations;
iter=0;

%Ratio between the last and the current objective terms
rapport=1;


%% Optimization loop

%Stopping criterion
while((rapport>epsilon))

    iter=iter+1;
    

    %Compute the temporal signals used to compute the gradient of the
    %regression term
    alpha_tonnext=ispecgram(spec_tonnext,fe,fen_ton,noverlap_ton);
    alpha_transnext=ispecgram(spec_transnext,fe,fen_trans,noverlap_trans);
    
    %Compute the gradients
    update_vector_ton=spectrogram(stepsize*(alpha_tonnext+alpha_transnext-f),fen_ton,noverlap_ton,taille_fenetre_ton,fe);
    update_vector_trans=spectrogram(stepsize*(alpha_tonnext+alpha_transnext-f),fen_trans,noverlap_trans,taille_fenetre_trans,fe);
   
    %Update the spectrogram thanks to the gradients
    spec_tonnext=update_opt(spec_tonnext,update_vector_ton);
    spec_transnext=update_opt(spec_transnext,update_vector_trans);
    

    
    
    
    %Shrinkage step    

    spec_trans = seuillage(spec_transnext,pond2vec*stepsize,p2,q2,group2);
    spec_ton = seuillage(spec_tonnext,pond1vec*stepsize,p1,q1,group1);


    
    fprintf('Iteration #%d\n',iter); 
    
    %Every 20 iterations, check if the objective have been majorly modified
    if((mod(iter,20)==0))
        
            %Save the last objective terms
            last_obj1 = obj1;
            last_obj2 = obj2;
            last_obj3 = obj3;

            %Compute the current objective terms
            alpha_ton=ispecgram(spec_ton,fe,fen_ton,noverlap_ton);
            alpha_trans=ispecgram(spec_trans,fe,fen_trans,noverlap_trans);
            [obj1 obj2 obj3] = compute_obj_function_signal(f,alpha_ton,alpha_trans,spec_ton,spec_trans,p1,q1,p2,q2,group1,group2);
            
            %Compute the sum of the ratio
            rapport=(abs(last_obj1-obj1)/last_obj1)+(abs(last_obj2-obj2)/last_obj2)+(abs(last_obj3-obj3)/last_obj3);
            
            fprintf('RATIO Iteration %d:  %f   \n',iter, rapport); 

    end;
    
    

    %Computation of parameters for the next iteration
    
    %Next computation points of the gradient
    spec_tonnext = compute_yk(t,spec_ton,spec_tonlast);
    spec_transnext = compute_yk(t,spec_trans,spec_translast);
    
    %Save the last computed values
    spec_tonlast = spec_ton;
    spec_translast = spec_trans;

    %Next fista t-parameter
    t = compute_tk(t);

    

end

%% Post process

%Compute the temporal signal for the three components from the spectrogram
alpha_ton=ispecgram(spec_ton,fe,fen_ton,noverlap_ton);
alpha_trans=ispecgram(spec_trans,fe,fen_trans,noverlap_trans);
bruit=f-alpha_ton -alpha_trans;


%Erase the 0 padding
alpha_ton=alpha_ton(1:long_sig);
alpha_trans=alpha_trans(1:long_sig);
bruit=bruit(1:long_sig);
f=f(1:long_sig);

%Draw the results
ih=spectrogram(alpha_ton, hanning(taille_fenetre_ton),floor(noverlap_ratio * taille_fenetre_ton), taille_fenetre_ton,1000);
it=spectrogram(alpha_trans, hanning(taille_fenetre_trans),floor(noverlap_ratio * taille_fenetre_trans), taille_fenetre_trans,1000);

Lsig=length(f);
axe_t = (1:Lsig)/fe;% Axe des temps
axe_ftrans= fe/2*linspace(0,1,taille_fenetre_trans/2);
axe_fton = fe/2*linspace(0,1,taille_fenetre_ton/2);

figure;
subplot(2,4,[1 4]);plot(axe_t,f);xlabel('Temps(sec)');ylabel('Amplitude(dB)');title('Resultat ');ylim([-1 1]);
subplot(2,4,5);plot(axe_t,alpha_trans);xlabel('Temps(sec)');ylabel('Amplitude(dB)');title('Transitoire estimée par l algorithme');ylim([-1 1]);
subplot(2,4,6);imagesc(axe_t,axe_ftrans,abs(it));xlabel('Temps(sec)');ylabel('Fréquence(Hz)');title('plan T-F représentant les transitoires');
subplot(2,4,7);plot(axe_t,alpha_ton);xlabel('Temps(sec)');ylabel('Amplitude(dB)');title('Tonale estimée par l algorithme');ylim([-1 1]);
subplot(2,4,8);imagesc(axe_t,axe_fton,abs(ih));xlabel('Temps(sec)');ylabel('Fréquence(Hz)');title('plan T-F représentant les tonales');


end
