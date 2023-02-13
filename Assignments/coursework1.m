%Generate two sinusoids x1(t), x3(t),  x5(t) with amplitude 2/p, 2/3p, 2/5p
% frequencies 300Hz and 500Hz. Add these three signals together and plot
% the results. Please explain the results.
% Note: x1 (t) has already been done in the class.
clear all;

%define the amplitudes and freqencies
a1 = 2/pi; a2 = 2/(3*pi); a3 = 2/(5*pi);
f1 = 100; f2 = 300; f3 = 500;

%define the sample timeframe
% f1 used to define period (T), this will display at least 4 cycles of x1
T = 1/f1;
Fs = 100000;
t = linspace(0,T*4, Fs);

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

%Discussion  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The waveform in Figure2 has multiple (local) minimums 
% and maximums, this is a result of adding multiple sinewaves together 
% that have different amplitudes and frequencies. 
% It can be seen from the Figure2 that waveform is beginning to take the
% form of a squarewave (e.g. we can see a half cycle where x is 0.4 to 0.6V
% between 0 and 0.005, with the other half at about -0.4 to -0.6V 
% from 0.005 to 1)

%Since a squarewave is an infinite series of sinewave harmonics, 
% we can approximate a squarewave by adding together multiple odd sinewaves
% with diminishing amplitudes and increasing frequencies.

%The fundamental frequency is 100Hz, this frequency will be a factor of
% each harmonic frequency. Since 2/3pi and 2/5pi are odd sinewaves (with
% 3*100 and 5*100 frequencies), these meet the requirement to approximate a
% square wave when combined, which is what we are seeing here. 
% Adding further harmonics (2/7pi(V) at 700Hz, 2/9pi(V) at 900Hz for
% example) will further steepen the rise/fall between each half-cycle, as
% well as flatten further the local minimums and maximums within each
% half-cycle – in other words, as we add more harmonics the wave will begin 
% to have a constant amplitude of 25V for 0.005s, and constant -25V for
% another 0.005s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add a constant dc value of 5V to the Sawtooth signal developed in class, 
% and explain the effect comparing to the original signal. 
% Compute the new mean, mean square and rms values.

%Code in the following block is taken from the lab, as stated in question
%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2

%Add a constant dc value of 5V to the Sawtooth signal developed in class, and
% explain the effect comparing to the original signal. 
% Compute the new mean, mean square and rms values.

% The code to generate the sawtooth was adapted from the lab code
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

% plot sawtooth wave
figure(1)
plot(t,v)
xlabel ('Time(s)')
ylabel ('Amplitude (V)')
title ('Figure 1: Sawtooth Signal')
ylim([0, 20])
grid on

%calculate the mean, mean square and root mean squared values
v_mean = mean(v);
v_mean_sq = mean(v.^2);
v_rms = sqrt(v_mean_sq);

% redefine waveform with DC offset
vp_off = ((E/T)*tp) + 5; % adding 5V for DC offset
v_off = [vp_off,vp_off,vp_off,vp_off]; % combine array vp four times to create 4 cycles

% plot DC-offset sawtooth wave
figure(2)
plot(t,v_off)
xlabel ('Time(s)')
ylabel ('Amplitude (V) with DC offset 5V')
title ('Figure 1: Sawtooth Signal')
ylim([0, 20])
grid on

%calculate the mean, mean square and root mean squared values with offset
v_off_mean = mean(v_off);
v_off_mean_sq = mean(v_off.^2);
v_off_rms = sqrt(v_off_mean_sq);

%{
 from figures (1) and (2) the effect of the 5V DC offset is visible,
 with every value v(t) shifted by 5V, hence the new maximimum point
 (amplitude) is 15V (previously it was 10V) and the new minimum point is
 5V (previously 0V)

 The DC component is also the mean of a waveform, hence adding 5V DC to
 this waveform, gives a new mean of v_mean+ 5 also

 mathematically, since the sawtooth wave is given as:
        v(t) = (E/T) * t [1]

 given that the mean of a waveform is defined as:
    v_mean = (1/T) * integral(v(t)^2) [2] (integrate w.r.t t, over 0:T)

 when [1] is subsituted [2] , we are left with:
    v_mean = E/2

 therefore adding a DC offset:
    v(t) = E0 + (E1/T) * t
    v_mean = E0 + (E1/2)

  hence the overall shift of 5V in this example
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% define the waveform for DC offset of -10V
v2 = zeros(1,arr_size) + a1 - 10; % create half cycle array with DC offset
vp2 = [v2, -v2]; % full cycle array
v2 = [vp2,vp2,vp2]; % array for three full cycles
t2 = linspace(0,T*3, length(v2)); % time sample of size v

%plot waveform
figure(1)
plot(t, v);
title('Squarewave');
xlabel('Time')
ylabel('Amplitude')
ylim([-30, 30])


%plot waveform
figure(2)
plot(t2, v2);
title('Squarewave with DC Offset (-10V)');
xlabel('Time')
ylabel('Amplitude')
ylim([-30, 30])


%calculate the mean, mean square and root mean squared values
v_mean = mean(v);
v_mean_sq = mean(v.^2);
v_rms = sqrt(v_mean_sq);

%calculate the mean, mean square and root mean squared values
v2_mean = mean(v2);
v2_mean_sq = mean(v2.^2);
v2_rms = sqrt(v2_mean_sq);


%{
This is as close to a mathematial 'proof' that matlab comments allow:

RMS for a waveform is given as below:
 rms = sqrt( (1/T) * integral(u(t)^2) ) [1]
 where:
   T is the period 
   u(t) is the wave function
for a bipolar square wave, we have:
    u1(t) = Vp for t = [0:t1]
    u2(t) = -Vp] for t = [t1:T]
        where [0:t1] is a half cycle

therefore, from eq.1 (integral wrt t)
    u1_rms = sqrt( (1/T) * integral(Vp^2) ) [2]
which reduces to:
    u1_rms^2  = (1/T) * (Vp^2) * (t1 - 0)
hence:
    u1_rms = Vp * sqrt(t1/T)
since t1/T is D, the duty cycle, and for a square wave with duty cycle
50%:
    u1_rms = Vp * sqrt(1) = Vp

Similarly, equation [2] can be adapted for u2_rms, to give:
    u2_rms^2  = (1/T) * ((-Vp)^2) * (T - t1)

with the rms for the full cycle then given as:
    u_rms = sqrt( u1_rms^2 + u2_rms^2 ) = Vp
%}

%redefine t and show mathemtically rms = E = 25
t = linspace(0,T, length(vp)); %redefine t for just one cycle
u = vp;
% from equations above, u_rms is given by:
u_rms = sqrt( (1/T) * trapz(t,u.^2) ); %trapz will approximate the integral

% verify with Matlab's rms
u_rms_mat = rms(u);
