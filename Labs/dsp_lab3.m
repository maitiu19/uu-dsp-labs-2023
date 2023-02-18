
%Qu2

% Qu3
Fs=100000;
Ts = 1 / Fs;
t = [0:1:250]*Ts;
f0=1000;
sig = 0.5;
for n=1:2:50
sig = sig+(2/(n*pi))*sin(2*pi*f0*n*t);
plot(t,sig)
pause(0.5)
end
xlabel('Time (seconds)')
ylabel('Signal Amplitude (V)')
 title (' Generating a Squarewave from its Fourier Components')
grid on 

% Qu5
Fs=10000;
t=[0:1:500]/Fs;
f0=100;
sig = pi/2;
for n=1:1:50
 sig=sig+(1/n)*(-1)^(n+1)*sin(2*pi*f0*n*t);
plot(t,sig)
pause(0.5)
end
xlabel('Time (seconds)')
ylabel('Signal Amplitude (V)')
title (' Generating a Sawtooth Signal its Fourier Components')
grid on


% Qu6
clear all
close all
Fs = 10000;
Ts = 1/Fs;
t=[0:1:399]*Ts;
%
for i=1:1:20
vp(i) = 10;
end
for i=21:1:100
vp(i)=0;
end
pulse = [vp,vp,vp,vp];
%
figure(1)
plot(t,pulse);
xlabel('Time(s)');
ylabel ('Amplitude(V)');
title('Periodic Pulse Signal')
grid on 

% Qu7
E = 10;
T = 10e-3;
tau = 2e-3;
DC = (E*tau)/T;
Fs=10000;
f0=100;
% no of samples per cycle = Fs/f0 = 10000/100 = 100
% 4 cycles => 4 x 100 = 400 samples
t=[0:1:399]/Fs;
sig= DC;
figure(2)
for n=1:1:50
 Cn = DC*exp(-j*pi*f0*n*tau)*sin(n*pi*tau/T)/(n*pi*tau/T);
 An = 2*abs(Cn);
 phin = angle(Cn);
 sig=sig +An*cos((2*pi*f0*n*t) + phin);
plot(t,sig)
pause (0.5)
end
xlabel('time (s)')
ylabel('Amplitude (V)')
title(' Synthesis of Periodic Pulse Signal by Adding its Fourier Components')
grid on

E = 10;
T = 10e-3;
tau = 2e-3;
f0 = 1/T;
spectrum(1) = E*tau/T; % set dc value
frequency(1) = 0;
for n = 2:1:50 % calculate harmonic values
 x = (n-1)*pi*tau/T;
 spectrum(n)= (E*tau/T)*(sin(x)/x)*exp(-j*n*pi*f0*tau);
 frequency(n)= (n-1)*f0;
end
figure(3)
stem (frequency, abs(spectrum))
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
title ('Amplitude Spectrum of Periodic Pulse Signal')
grid on
figure(4)
stem (frequency, angle(spectrum))
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
title ('Phase Spectrum of Periodic Pulse Signal');
grid on