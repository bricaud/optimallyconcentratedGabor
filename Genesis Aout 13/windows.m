% do the examples from Guillaume
%


load window;
load window_1;
load window_2;

% first window
%

figure;
sgram(window.gamma,'lin','nf','tc');

[Lpad1,a1,M1,lt1]=optimalsampling_new(window.Len,2,window.d,window.s)
b1=Lpad1/M1


% second window
%

figure;
sgram(window1.gamma,'lin','nf','tc');

[Lpad2,a2,M2,lt2]=optimalsampling_new(window1.Len,2,window1.d,window1.s)
b2=Lpad2/M2



% third window
%

figure;
sgram(window2.gamma,'lin','nf','tc');

[Lpad3,a3,M3,lt3]=optimalsampling_new(window2.Len,2,window2.d,window2.s)
b3=Lpad3/M3
