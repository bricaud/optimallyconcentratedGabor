% more questions about the speed of the LTFAT
%

L=367500;
f=rand(L,1);
g=rand(L,1);

a=1;
M=30;

% slow computation on time side
tic; dgt(g,f,a,M,'lt', [1 2]); toc

% fast computation on frequency side
tic; dgt(g,f,L/M,L/a, 'lt', [1 2]);toc

% question: why is there the big difference between computation?
% is this connected to what you call parameters c and d?
