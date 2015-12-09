function imagescMatrixOnNonSeparableLattice(matrix,Lt,Tmax,Fmax)

%imagescMatrisOnNonSeparableLattice 
%       Plot matrix computed on a non separable lattice (for example a non
%       separable DGT)
%
%       Inputs:
%           Matrix: matrix to plot
%           Lt: parameters of the lattice, in the same format than for the
%              DGT computation (see MATRIX2LATTICETYPE for more details)
%           Tmax: last label on the x axis
%           Fmax: last label on the y axis for the first column of the
%              matrix 
M = size(matrix,1);
L = size(matrix,2);
Lt1 = Lt(1);
Lt2 = Lt(2);

% Prepare data : repeat data in order to have "longer" (and overlaping)
% frequency bins
matrixTemp = zeros(Lt2*M,L);

for i = 1:M
    matrixTemp((i-1)*Lt2+1:i*Lt2,:) = repmat(matrix(i,:),Lt2,1);
end

% Shift data to overlap bins and get final matrix to be plotted
matrixToPlot = NaN(Lt2*(M+1),L);

M2 = Lt2*(M+1);
for i = 1:L,
    m = mod(Lt1*(i-1),Lt2);
    matrixToPlot(m+1 : M2-(Lt2-m),i) = matrixTemp(:,i);
end
% Faire les axes
axisX = 0:Tmax/(L-1):Tmax;
axisY = 0:Fmax/((M-1)*Lt2):Fmax+Fmax/((M-1)*Lt2);

set(0,'DefaultAxesFontSize',18)
imagesc(axisX, axisY, matrixToPlot)
axis xy;
m = max(max(matrix)); caxis([m-60 m]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
