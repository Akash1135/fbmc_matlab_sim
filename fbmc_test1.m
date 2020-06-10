clc
close all
clear all
N = input("Enter the value")


H1=0.971960;

H2=sqrt(2)/2;

H3=0.235147;

factech=1+2*(H1+H2+H3);

hef(1:4*N)=0;

for i=1:4*N-1

   hef(1+i)=1-2*H1*cos(pi*i/(2*N))+2*H2*cos(pi*i/N)-2*H3*cos(pi*i*3/(2*N));

end

hef=hef/factech;
% Prototype filter impulse response

h=hef;  

 

% Initialization for transmission

Frame=1;

y=zeros(1,4*N+(Frame-1)*N/2);

 

s=zeros(N,Frame);

for ntrame=1:Frame
    % OQAM Modulator

    if rem(ntrame,2)==1

        s(1:2:N,ntrame)=sign(randn(N/2,1));

        s(2:2:N,ntrame)=j*sign(randn(N/2,1));

    else

        s(1:2:N,ntrame)=j*sign(randn(N/2,1));

        s(2:2:N,ntrame)=sign(randn(N/2,1));

    end
    x=ifft(s(:,ntrame));

 

% Duplication of the signal

x4=[x.' x.' x.' x.'];

 

% We apply the filter on the duplicated signal

signal=x4.*h;

%signal=x4;

% Transmitted signal

y(1+(ntrame-1)*N/2:(ntrame-1)*N/2+4*N)=y(1+(ntrame-1)*N/2:(ntrame-1)*N/2+4*N)+signal;

 

end
yfft = [zeros(1,length(y)/2) fft(y) zeros(1,length(y)/2)]

yTx = ifft(yfft);

yTxFreq = fft(yTx,8192);

yTxFreqAbs = abs(yTxFreq);

yTxFreqAbs = yTxFreqAbs/max(yTxFreqAbs);

yTxFreqAbsPwr = 20*log(yTxFreqAbs);

 

figure(); plot(h);xlim([0 length(h)]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('Amplitude')
title('1.Co-efficient of prototype filter')

figure(); plot(s);xlim([-1.2 1.2]);ylim([-1.2 1.2]);set(gca,'yticklabel',[]);
xlabel('In-Phase ')
ylabel('Quadrature ')
title('2.Constellation Diagram')

figure(); stem(real(s));xlim([1 length(s)]);ylim([-1.2 1.2]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('Amplitude')
title('3.In-Phase data to be transmitted')

figure(); stem(imag(s));xlim([1 length(s)]);ylim([-1.2 1.2]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('Amplitude')
title('4.Quadrature data to be transmitted')

figure(); plot(real(x));xlim([1 length(x)]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('IFFT(plot(3))')
title('5.IFFT of in-phase components')

figure(); plot(imag(x));xlim([1 length(x)]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('IFFT(plot(4))')
title('6.IFFT of Quadrature components')

figure(); plot(real(x4));xlim([0 length(x4)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
xlabel('time')
ylabel('4 * Amplitude ')
title('7.Duplication of in-phase component')

figure(); plot(imag(x4));xlim([0 length(x4)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
title('8.Duplication of Quadrature component')
xlabel('time')
ylabel('4 * Amplitude ')

figure(); plot(real(signal));xlim([0 length(signal)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
title('9.Prototype Filter * Duplication of inphase Component')
xlabel('time')
ylabel('Amplitude')

figure(); plot(imag(signal));xlim([0 length(signal)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
title('10.Prototype Filter * Duplication of Quadrature Component')
xlabel('time')
ylabel('Amplitude')

figure(); plot(real(y));xlim([0 length(y)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
title('11.Frame 1 similar to 9') 
xlabel('time')
ylabel('Amplitude')

figure(); plot(imag(y));xlim([0 length(y)]);ylim([-0.5 0.5]);set(gca,'yticklabel',[]);
title('12.Frame 1 similar to 10') 
xlabel('time')
ylabel('Amplitude')

figure(); plot(real(yfft));xlim([0 length(yfft)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('13.FFT of 12: inphase * h(t)')
xlabel('time')
ylabel('Amplitude')

figure(); plot(imag(yfft));xlim([0 length(yfft)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('14.FFT of 13: Quadrature * h(t)')
xlabel('time')
ylabel('Amplitude')

figure(); plot(abs(yfft)/max(abs(yfft)));xlim([0 length(yfft)]);set(gca,'yticklabel',[]);%ylim([0 1]);
title('15.FFT result with zero padding 11 and 12') 
xlabel('Frequency')
ylabel('Magnitude')

figure(); plot(20*log(abs(yfft)/max(abs(yfft))));xlim([0 length(yfft)]);set(gca,'yticklabel',[]);%ylim([-10 0]);
title('16.Magnified view of passband') 
xlabel('time')
ylabel('Amplitude')

figure(); plot(real(yTx));xlim([0 length(yTx)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('17.Real number of IFFT result')
xlabel('time')
ylabel('Amplitude')

figure(); plot(imag(yTx));xlim([0 length(yTx)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('18.Imaginary number of IFFT result')
xlabel('time')
ylabel('Amplitude')

figure(); plot(abs(yTx)/max(abs(yTx)));xlim([0 length(yTx)]);ylim([0 1]);set(gca,'yticklabel',[]);
title('19.Magnitude of IFFT result of 15:FFT with zero padding')
xlabel('time')
ylabel('Amplitude')

figure(); plot(20*log(abs(yTx)/max(abs(yTx))));xlim([0 length(yTx)]);ylim([-100 0]);set(gca,'yticklabel',[]);
title('20.Representation of IFFT result in db')
xlabel('time')
ylabel('Amplitude')

figure(); plot(real(yTxFreq));xlim([0 length(yTxFreq)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('21.Real number after FFT with zero padding')
xlabel('time')
ylabel('Amplitude')

figure(); plot(imag(yTxFreq));xlim([0 length(yTxFreq)]);set(gca,'yticklabel',[]);%ylim([-0.5 0.5]);
title('22.Imaginary number after FFT with zero padding')
xlabel('time')
ylabel('Amplitude')
figure(); plot(yTxFreqAbs);xlim([0 length(yTxFreqAbs)]);ylim([0 1]);set(gca,'yticklabel',[]);
title('23.Magnitude of real number & imaginary ')
xlabel('time')
ylabel('Magnitude')
figure();plot(yTxFreqAbsPwr); xlim([0 length(yTxFreqAbsPwr)]);ylim([-100 0]);set(gca,'yticklabel',[]);
title('24.Magnitude of real & imaginary number in decibel(db)')
xlabel('time')
ylabel('Magnitude')