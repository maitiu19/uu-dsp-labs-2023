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

clear all;
close all;
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
Question 3
***********************************************************************
***********************************************************************
%}

clear all;

an = [1, 2, 4, -9, 1];
bn = [2, -1, 3, 3, 0];

cn = 2 * an .* bn;
dn = an + bn;
en = an/2;



%{
Question 4
***********************************************************************
Pulse wave generation code below was adapted from the lab3 file and
modified as required, Fourier part coded separately below
***********************************************************************
%}

clear all;
close all;

T = 2*pi;   % signal period
t= linspace(-T, T, 400);  % 400 samples between -2pi and 2pi

% define the amplitude array for one full cycle, T = 2pi
vp = zeros(1,200);  % populate with zeros
vp(101:200) = 1;    % second half cycle
pulse = [vp,vp];    % define v for full sample range, 4pi 

figure(1)
plot(t,pulse);
xlabel('Time(s)');
ylabel ('Amplitude(V)');
title('Periodic Pulse Signal from -2\pi to 2\pi');
set(gca,'XTick',-2*pi:pi:2*pi);
ax.XTickLabel = {'-2\pi','-\pi','0','\pi','2\pi'};
ylim([-.5, 1.5]);
xlim([-2*pi, 2*pi]);
grid on

%{
For the Fourier part:

This is working as a visual only, there is clearly some fundamental error
due in the code below, or in the math or error or both

Since this is an odd funciton, A_n will be zero, so was not included, but
when it was included, it evaluated to a non-zero value, again, there is an
error in here somewhere, it is of magnitude x10^(-15) so very small, could
be due to rounding error and trapz method inaccuracy

It might be related to the fact that I am taking T to be 2pi, but
evaluating over a period of 4pi (-2pi to +2pi) - adding 2*T to the
equations for b_n and sig have more issues are N gets large

2*T I think is justified mathematically (?), but that would mean 4T for A0
which gives a value of 0.025, when (mean(vp) = 0.5, so I left this as T

The main issue when N is small, <10, it is just a sine wave,
when it is >10 it begins to resemble the function, 10
components in this sequence should have a stronger resemblance to the true
fucntion than a simle sine wave! Going larger, e.g. N>200, there amplitude
is distorted, perhaps signal noise, but no idea where these errors are
caused initially.

trapz should be acureate enough (?) not to result in a completely different
waveform at N= 1000, so not sure where the error is
%}

clear all;
close all;

E = 1;
T = 2*pi;   % signal period
t= linspace(-T, T, 400);  % 400 samples between -2pi and 2pi
N = 200;
% define the amplitude array for one full cycle, T = 2pi
vp = zeros(1,200);  % populate with zeros, f(0:pi) is off, i.e. 0
vp(101:200) = E;    % first half cycle is 'on', e.g. f(-pi:0) = 1
vp = [vp, vp];

% Evaluate the Fourier series over the full period
a0 = 1/(2*T) * trapz(t, vp); % should be T*4 ??
an = zeros(1,N);
bn = zeros(1,N);
for n = 1:N
    an(n) = 1/(T) * trapz(t, (vp .* cos((n * pi .* t)/(T))));
    bn(n) = 1/(T) * trapz(t, (vp .* sin((n * pi .* t)/(T))));
end


% Setting up the Fourier derived signal, a_n is removed since it should
% evaluate to zero anyway.
sig = zeros(size(t));
sig = sig + a0;
for n = 1:N
    sig = sig  + bn(n) * sin((n * pi .* t)/(T));
end

% Plot the original signal and the reconstructed signal
figure(2);
plot(t, vp);
hold on;
plot(t, sig, '--');
title('True Signal vs Fourier Dervived Signal');
legend('True Signal', 'Fourier Dervived Signal');
xlabel('Time(s)');
ylabel ('Amplitude(V)');
set(gca,'XTick',-2*pi:pi:2*pi);
ax.XTickLabel = {'-2\pi','-\pi','0','\pi','2\pi'};
xlim([-2*pi, 2*pi]);
grid on


%{
************************************************************************
The following block of code was taken from the week 3 lab, Qu.7, but
adapted for the function described in this question

I've included this since this is more reliable than the method used above
which has some unknown errors somewhere, plus I realised after the week 4
lab maybe the requirement was to use Complex Fourier for this question
************************************************************************
%}

clear all;
close all;

E = 1;          % amplitude is 1V
T = 2*pi;       % period cycle from -pi to pi, therefore T = 2pi
tau = pi;       % tau/T is the on/off ratio, this pulse is 50/50
DC = (E*tau)/T; % Duty Cycle
Fs= 2*T;        % sampling over two periods
f0= 1/T;

% 500 sample points over two cycles, -2pi to 2pi
t = linspace(-T, T, 500);
sig= DC;

% define the function's  Fourier series components and plot as each
% component is added to the signal
figure(2)
for n=1:1:50
    Cn = DC*exp(-j*pi*f0*n*tau)*sin(n*pi*tau/T)/(n*pi*tau/T);
    An = 2*abs(Cn);
    phin = angle(Cn);
    sig=sig +An*cos((2*pi*f0*n*t) - phin);
    plot(t,sig)
    pause (0.5)
end
xlabel('time (s)')
ylabel('Amplitude (V)')
title(' Synthesis of Periodic Pulse Signal by Adding its Fourier Components')
grid on