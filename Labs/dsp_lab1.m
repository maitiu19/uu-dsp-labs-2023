A = [1,2,3; 4,5,6;,7,8,9];
B = [10, 11, 12, 13; 14, 15, 16, 17; 18, 19, 20, 21];
C = [9, 8, 7; 6, 5, 4; 3, 2, 1];

% A squared
A_sq = A.^2;

%A transposed
A_transp = A.';

%Solve the following linear equations by converting to matrices
% 3x1 + x2 = 18; 4x1 + 2x2 = 21

%transform into an AX = B form
A = [3, 1; 4, 2];
B = [18; 21];

%use linsolve to find x
X = linsolve(A, B);

%solve for x using equtomatrix
syms x1 x2 x3
eq1 = 5*x1 - 3*x2 - 2*x3 == 31;
eq2 = 2*x1 + 6*x2 + 3*x3 == 4;
eq3 = 4*x1 + 2*x2 - x3 == 30;

[A, B] = equationsToMatrix([eq1, eq2, eq3], [x1, x2, x3]);

X = linsolve(A, B);

% create an array "t" of 100 data points from 0, increasing linearly with
% an increment of 1. Print or plot the array "t" to the screen
t = linspace(0,99);

%a discrete sinusoid consisting of 100 data points with 
% amplitude 2/π, 
% sample period T = 1ms,
% frequency 100Hz
x = linspace(0,1,100);
x1 = (2/pi) * (sin(2*pi*100*x*1));

% Generate two additional sinusoids x3(t) and x5(t) with amplitude 2/3π 
% and 2/5π and frequencies 300 Hz and 500Hz. Add x1(t), x3(t) and x5(t) 
% together and plot the result.
x2 = (2/3*pi) * (sin(2*pi*300*x*1));
x3 = (2/5*pi) * (sin(2*pi*500*x*1));

%plot all three waves on the one figure
figure
plot(x, x1, x, x2, x, x3);
title('Sine waves of first 3 harmonics');
xlabel('Time')
ylabel('Amplitude')