%{
Question 1
***********************************************************************
A sinusoidal signal of amplitude 10V and frequency 1kHz is sampled at a
frequency f s = 10kHz, beginning at time t = 0. Determine the first ten
values of the sampled sequence.

Since there is a sample taken once every 0.1 milliseconds (10kHz) and the 
signal cycle completes once every 1 millisecond (1kHz), the number of
samples taken will be 10, over a 10 millisecond period starting from zero

***********************************************************************
%}

clear all
close all
fs = 10000;                 % 10kHz, sample frequency
f0 = 1000;                  % 1kHz, signal frequency
E = 10;                     % 10V, amplitude
T = 1/f0;                   % signal time period
Ts = 1/fs;                  % Ts = sample interval ( period )

tp = 0:Ts:T;              % array of sample times for one period
sample_vals = E * sin(2*pi*f0*tp); % get the sample values, i.e. x(t)


t = linspace(0, T, fs);     % set up the time array for plot
sig = E * sin(2*pi*f0*t);   % signal array for plot

% plot for clarity/sanity check
figure(1)
plot(t,sig, tp, sample_vals, 'o', 'MarkerFaceColor', 'b');
ylim([-11, 11]);

% output the values for the first 10 samples in the sequence
disp(sample_vals);



%{
Question 2
***********************************************************************
Show that a sinusoid of amplitude 10V and frequency 11kHz sampled at
f s = 10kHz as in Q.1 above is an alias of the 1kHz sampled signal.
***********************************************************************
%}

alias_fs = 10000;                 % 10kHz, sample frequency
alias_f0 = 11000;                  % 11kHz, signal frequency
alias_E = 10;                     % 10V, amplitude
alias_T = 1/alias_f0;                   % signal time period
alias_Ts = 1/alias_fs;                  % Ts = sample interval ( period )
alias_tp = 0:alias_Ts:alias_T*11;              % array of sample times for one period
alias_sample_vals = E * sin(2*pi*f0*alias_tp); % get the sample values, i.e. x(t)
alias_t = linspace(0, alias_T*11, alias_fs);     % set up the time array for plot
alias_sig = E * sin(2*pi*alias_f0*alias_t);   % signal array for plot

% plot for clarity/sanity check
fig2 = figure(2);
plot(t,sig);
hold('on');
plot(alias_t, alias_sig, 'g-.')
plot(tp, sample_vals, 'kx', 'MarkerSize',16)
plot(tp, alias_sample_vals,'d','MarkerFaceColor', 'k');
ylim([-11, 11]);
xlabel('Time (seconds)')
ylabel('Signal Amplitude (V)')
title ('11kHz Alias Signal vs. 1kHz Signal')
legend('1kHz Signal', '11kHz Alias Signal', '1kHz Sample Points', '11kHz Alias Sample Points')


% output the values for the first 10 samples in the alias sequence
disp(alias_sample_vals);


%{
Question 4
***********************************************************************
Pulse wave generation code below was adapted from the lab3 file and
modified as required
***********************************************************************
%}

clear all;
close all;

T = 2*pi;   % signal period
t= linspace(-T, T, 400);  % 400 samples between -2pi and 2pi

% define the amplitude array for one full cycle, T = 2pi
vp = zeros(1,200);  % populate with zeros, f(-pi:0) is off, i.e. 0
vp(101:200) = 1;    % second half cycle is 'on', e.g. f(0:pi) = 1
pulse = [vp,vp];    % define v for full sample range, 4pi

figure(1)
plot(t,pulse);
xlabel('Time(s)');
ylabel ('Amplitude(V)');
title('Periodic Pulse Signal');
set(gca,'XTick',-2*pi:pi:2*pi);
ax.XTickLabel = {'-2\pi','-\pi','0','\pi','2\pi'};
ylim([-.2, 1.2]);
xlim([-2*pi, 2*pi]);

clear all;

% Qu7
E = 1;
T = 2*pi;
tau = pi;
DC = (E*tau)/T;
f0=1/T;
% no of samples per cycle = Fs/f0 = 10000/100 = 100
% 4 cycles => 4 x 100 = 400 samples
t= linspace(-T, T, 400);  % 400 samples between -2pi and 2pi
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
set(gca,'XTick',-2*pi:pi:2*pi);
ax.XTickLabel = {'-2\pi','-\pi','0','\pi','2\pi'};
ylim([-.2, 1.2]);
xlim([-2*pi, 2*pi]);



