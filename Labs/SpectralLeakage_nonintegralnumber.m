fs = 1000;
T = 1/fs;
A=0.5;
f=2.5;
N = 400; % Create a signal with N samples = 1 second
n = 0:N-1; % sample numbers
t = n*T; % sampling times over a one second duration.
x = A*sin(2*pi*6*t);
x = [x zeros(1,3000)];
N = length (x);
n = 0:N-1;


 
%plot the magnitude spectrum 
X_mags = abs(fft(x));
 
subplot(2,1,1)
plot(x)
xlabel('Samples');
ylabel('Amplitude')
title('Time-Domain Signal');
 
num_disp_bins = 10;
subplot(2,1,2)
stem([0:num_disp_bins-1], X_mags(1:num_disp_bins));
hold on
stem([0:num_disp_bins-1], X_mags(1:num_disp_bins),'k.');
hold off
xlabel('Frequency Bins');
ylabel('Magnitude');
title('Frequency Content Magnitudes');



%Illustrate Basis functions

plot(n,x,'k')
xlabel('Samples');
ylabel('Amplitude')
for k = 0 : 20
    plot(x,'kx')  
    hold on
    cos_basis = A*cos(2*pi*k*n/N);
    sin_basis = A*sin(2*pi*k*n/N);
     
    plot(n,cos_basis,'r')
    plot(n,sin_basis,'g')
    title({['Signal being anlaysed (black) with basis functions'...
        ' that have ' num2str(k) ' cyles over '  ...
        num2str(N) ' samples']...
        ['Cosine basis function shown in red, Sine in green']})
    hold off
    pause
end