
<<<<<<< HEAD
%Qu1
clear all 

% Fourier Series

Fs = 100000; % Sampling frequency = 100 kHz
t = [0:1:250]/Fs; % Constructs a time array t with 251 samples
% increasing from 0 to 250/100000 (2.5ms).

y1=(2/pi)*sin(2000*pi*t) ; % creates array y1 of 251 samples of 1kHz
% fundamental component (1st harmonic)
plot (t,y1);
xlabel('time (seconds)') 
ylabel('signal amplitude')

%plots fundamental component

% creates array y2 of 251 samples of
% 3kHz 3rd harmonic component
y2=(2/(3*pi))*sin(6000*pi*t); 
plot(t,y2) % plots 3rd harmonic component

y=y1+y2; % adds 1st and 3rd harmonic components
plot(t,y) % plots 1st and 3rd harmonic components

% creates array y3 of 251 samples of
% 5kHz 5th harmonic component
y3=(2/(5*pi))*sin((10000*pi*t) + pi/2) ;
plot(t,y3)

y=y1+y2+y3; % adds 1st , 3rd and 5th harmonic components
plot(t,y) % plots 1st , 3rd and 5th harmonic components

% creates array y4 of 251 samples of
% 7kHz 7th harmonic component
y4=(2/(7*pi))*sin(14000*pi*t); 

y=y1+y2+y3+y4; % adds first, 3rd, 5th and 7th harmonic components
plot(t,y) % plots 1st , 3rd, 5th and 7th harmonic components

% creates array y5 of 251 samples of
% 9kHz 9th harmonic component
y5=(2/(9*pi))*sin(18000*pi*t);

% creates array y6 of 251 samples of
% 11kHz 11th harmonic component
y6=(2/(11*pi))*sin(22000*pi*t);

% adds 1st, 3rd, 5th, 7th, 9th and 11th harmonic components
y=y1+y2+y3+y4+y5+y6;
plot(t,y) % plots 1st , 3rd, 5th , 7th , 9th and 11th harmonic components

% creates array y7 of 251 samples of
% 13kHz 13th harmonic component
y7=(2/(13*pi))*sin(26000*pi*t);
plot(t,y7)

% creates array y8 of 251 samples of
% 15kHz 15th harmonic component
y8=(2/(15*pi))*sin(30000*pi*t);

% creates array y9 of 251 samples of %17kHz 17th harmonic component
y9=(2/(17*pi))*sin(34000*pi*t);

% creates array y10 of 251 %samples of 19kHz 19th harmonic component
y10=(2/(19*pi))*sin(38000*pi*t);

% adds 1st , 3rd , 5th , 7th , % 9th , 11th , 13th , % 15th , 17th,
% and 19th harmonic components
y=y1+y2+y3+y4+y5+y6+y7+y8+y9+y10;
plot(t,y) 

xlabel('time (seconds)') % labels x-axis of amplitude-time graph
ylabel('signal amplitude') % labels y-axis of amplitude-time graph
grid on

% adds 1st , 3rd , 5th , 7th , % 9th , 11th , 13th , % 15th , 17th,
% and 19th harmonic components
y= 0.5 + y1+y2+y3+y4+y5+y6+y7+y8+y9+y10;
plot(t,y) 

xlabel('time (seconds)') % labels x-axis of amplitude-time graph
ylabel('signal amplitude') % labels y-axis of amplitude-time graph
grid on


%Qu3
clear all
=======
%Qu2

% Qu3
>>>>>>> 0e33a23b5d29071b92dd749e54775320383b4242
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


% Qu4
clear all

Fs=10000;
t=[0:1:500]/Fs;
f0=100;
sig = 0;
for n=1:1:50
sig=sig+(1/n)*(-1)^(n+1)*sin(2*pi*f0*n*t);
plot(t,sig)
pause(0.5)
end
xlabel('Time (seconds)')
ylabel('Signal Amplitude (V)')
title (' Generating a Sawtooth Signal its Fourier Components')
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
vp(i) = 1;
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