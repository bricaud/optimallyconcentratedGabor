function [ groups ] = create_groups_bands( g_l, mat_h, mat_l, fs )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [ groups ] = create_groups_bands( g_l, mat_h, mat_l, fs )
%   
% 
% FUNCTION:
%   Create groups in a spectrogram in order to be used by mixed norms 
% 
% 
% INPUT:
%       g_l: Width of a group
%       mat_h, mat_l: Size of the spectrogram
%       fs: Sample frequency
%
% OUTPUT:
%       groups   : Matrix of the same size than the spectrogram, which
%       contains for each coefficients, the group to which it belongs
%       
%
% NOTA: Frequency bands width are arbitrary defined (1100 Hz) but should be of
% equal for all bands (in order to not introduce groups of different size for regularization) 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define bands
bands=0:1100:fs/2+1100;
tab_entier=[1:size(bands,2)];

%Number of groups by line
gr_by_line=ceil(mat_l/g_l);
%Number of groups by column
gr_by_col=size(bands,2)-1;
%Number of groups
group_nb=gr_by_line*gr_by_col;

groups=cell(group_nb,1);

%For each coefficient, computing to which group it belongs
for iter = 1 : mat_h*mat_l,
    num_col=floor((iter-1)/mat_h)+1;
    num_line=iter-mat_h*(num_col-1);
    freq_line=(num_line)*fs/(2*mat_h);
    num_bands=max(tab_entier(bands<freq_line));
    num_gr=(num_bands+gr_by_col*floor((num_col-1)/g_l));
    groups{num_gr}=[groups{num_gr} iter];
end
end
