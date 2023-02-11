


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a discrete sinusoid consisting of 100 data points with 
% amplitude 2/π, 
% sample period T = 1ms,
% frequency 100Hz
x1 = (2/pi) * (sin(2*pi*100*t*1));

% Generate two additional sinusoids x3(t) and x5(t) with amplitude 2/3π 
% and 2/5π and frequencies 300 Hz and 500Hz. Add x1(t), x3(t) and x5(t) 
% together and plot the result.
x2 = (2/(3*pi)) * (sin(2*pi*300*t*1));
x3 = (2/(5*pi)) * (sin(2*pi*500*t*1));

%plot all three waves on the one figure
figure(1)
plot(t, x1, t, x2, t, x3);
title('Sine waves of first 3 harmonics');
xlabel('Time')
ylabel('Amplitude')

%plot all three waves on the one figure
x = x1 + x2 + x3;
figure(2)
plot(t, x);
title('Sine waves of first 3 harmonics');
xlabel('Time')
ylabel('Amplitude')