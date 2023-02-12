%Generate two sinusoids x1 (t) , x3 (t) and x5 (t) with amplitude 2/p, 2/3p and 2/5p
% and frequencies 300Hz and 500Hz. Add these three signals together and plot
% the results. Please explain the results.
% Note: x1 (t) has already been done in the class.

clear all;

%define the amplitudes and freqencies
a1 = 2/pi; a2 = 2/(3*pi); a3 = 2/(5*pi);
f1 = 100; f2 = 300; f3 = 500;

%define the sample timeframe
% f1 is used to define period (T) to display at least 4 cycles of x1
T = 1/f1;
n = 100000;
Ts = 1/n;
t = linspace(0,T*4, n);

% define the waveforms
x1 = a1 * (sin(2*pi * f1 * t));
x2 = a2 * (sin(2*pi * f2 * t));
x3 = a3 * (sin(2*pi * f3 * t));


%plot all three waves on the one figure
figure(1)
plot(t, x1, t, x2, t, x3);
title('Sine waves with different Amplitude and Freq.');
xlabel('Time')
ylabel('Amplitude')

% add the waves together to generate combined waveform, x
x = x1 + x2 + x3;

% plot the combined waveform
figure(2)
plot(t, x);
title('Combined Waveform');
xlabel('Time')
ylabel('Amplitude')


%Add a constant dc value of 5V to the Sawtooth signal developed in class, 
% and explain the effect comparing to the original signal. 
% Compute the new mean, mean square and rms values.

%Code in the following block is taken from the lab, as stated in question
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sawtooth
clear all
E = 10; % E = signal amplitude
T = 0.01; % period of signal = T; f = 1/T
Fs = 100000; % Fs = sample frequency
Ts = 1/Fs; % Ts = sample interval ( period )
tp = [0:Ts:T]; % generate time array of duration one period
vp = (E/T)*tp; % generate one cycle of the sawtooth signal in array vp
v = [vp,vp,vp,vp]; % combine array vp four times to create 4 cycles
t=[0:1:(4*T/Ts)+4]*Ts; % ensure that length of time array t equates with that of voltage array v
figure(1)
plot(t,v)
xlabel ('Time(s)')
ylabel ('Amplitude (V)')
title ('Figure 1: Sawtooth Signal')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%

vp = ((E/T)*tp) + 5; % adding 5V for DC offset
v = [vp,vp,vp,vp]; % combine array vp four times to create 4 cycles
t=[0:1:(4*T/Ts)+4]*Ts; % ensure that length of time array t equates with that of voltage array v
figure(2)
plot(t,v)
xlabel ('Time(s)')
ylabel ('Amplitude (V) with DC offset 5V')
title ('Figure 1: Sawtooth Signal')
grid on


%Sqaure wave
clear all;
%define the amplitudes and freqencies
a1 = 25;
f1 = 100;

%define the sample timeframe
% f1 is used to define period (T) to display at least 4 cycles of x1
T = 1/f1;
Fs = 100000;
Ts = 1/Fs;

% define the waveform
arr_size = int64((T/(2*Ts))); %size of array for half cycle
v1 = zeros(1,arr_size) + a1; % create half cycle array with amplitude 25
vp = [v1, -v1]; % full cycle array
v = [vp,vp,vp]; % array for three full cycles
t = linspace(0,T*3, length(v)); % time sample of size v

%plot waveform
figure(1)
plot(t, v);
title('Squarewave');
xlabel('Time')
ylabel('Amplitude')
ylim([-30, 30])


% define the waveform for DC offset of -10V
v1 = zeros(1,arr_size) + a1 - 10; % create half cycle array with DC offset
vp = [v1, -v1]; % full cycle array
v = [vp,vp,vp]; % array for three full cycles
t = linspace(0,T*3, length(v)); % time sample of size v


%plot waveform
figure(2)
plot(t, v);
title('Squarewave with DC Offset (-10V)');
xlabel('Time')
ylabel('Amplitude')
ylim([-30, 30])