fs=1000;
t=0:1/fs:1-1/fs;
f1=2.5;
x1=0.5*cos(2*pi*f1*t +0.2); 
figure; plot(x1,'r'); hold on
% x1 = x1.*hanning(length(x1))'; %Length of Hanning 
                              %window needs to be same as x1
figure 
hanwin=hanning(1000); %just to visualise the 
                     %shape of Hanning window
plot(hanwin,'b')
hold on
plot (x1, 'g') %after x1 has been multiplied 
                %by Hanning window
plot (0.5*cos(2*pi*f1*t +0.2), 'g')
% x1 = [x1 zeros(1, 11000)];  
subplot(2,1,1)
plot(x1)
xlabel('samples');
ylabel('amplitude');
X1=fft(x1);
subplot(2,1,2)
plot([0:length(X1)-1], abs(X1))
ylabel('Magnitude')
xlabel('Bins')


