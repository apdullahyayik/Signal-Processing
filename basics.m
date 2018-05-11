% Some useful basic signal processing techniques without using SP tool in Matlab.
%DOT PRODUCT
%TIME DOMAIN CONVULUTION
%DISCRETE FOURIER TRANSFORM
%SHORT-TIME FOURIER TRANSFORM
%TAPERING AND MULTI-TAPERING
%FREQUENCY DOMAIN CONVOLUTION
%CONVOLITION AS A FILTER
%MORLET WAVELET
%COMPLEX WAVELET
%WAVELET CONVOLUTION
%
%
%
%                                  Authored by , Apdullah Yayık 2017


%% Dot product
len=1000;
frequency=50;
time = (0: len-1)/len; 

signal=rand(1, len);
figure, plot(signal)

sine_wave =2*sin(2*pi* 10 *time + 0);
figure, plot(sine_wave)

signaldotkernel=signal.*sine_wave;
figure, plot(signaldotkernel);

%% Sum of Sine Waves

len=1000;
time = (0: len-1)/len; 
sinwave1=2*sin(2*pi* 10 *time + 0);
sinwave2=2*sin(2*pi* 20 *time + 0);
sumofsinwaves =2*sin(2*pi* 10 *time + 0)+2*sin(2*pi* 20 *time + 0);
figure, plot(sinwave1, 'b'), hold on, plot(sinwave2, ''), hold on,  plot(sumofsinwaves,'r'), 

%% Convolution Time Domain
% sinyal üzerinde kernel i gezdirerek dot product yap

load('inputData.mat')
signal=inputData(1,:); % 1 electrode EEG data
%  zeroadd=zeros(1, 50);
%  signazeroadded=[signal zeroadd];
kernel1=(50:-1:1)./50;
kernel2=hann(50)';
c=zeros(1, 2500); 
for i=1:2500
c(i)=sum(signal(i:i+49).*kernel2);
end

figure, plot(signal)
hold on, plot(c, 'r-')

 % Build-in Function
close all,
c2=conv(signal, kernel2);

figure, plot(signal)
hold on, plot(c, 'r-')

%% Convolution Freq Domain
% fourier of signal * fourier of kernel (point by point)

% Fourier of signal 
close all; clear all; clc
load('inputData.mat')
data = inputData(1,:); 
N = length(data); 
Fs=1000;
f = Fs/2*linspace(0,1,N/2+1);

% initialize Fourier coefficients
fourierofSignal = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform
for fi=1:N
% create sine wave

% sine_wave =2*sin(2*pi* (fi-1) *time + 0);
sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
fourierofSignal(fi) = sum(sine_wave.*data);
end

% Fourier of kernel
data = hann(50)';
N = length(data); 
Fs=1000;
f = Fs/2*linspace(0,1,N/2+1);

% initialize Fourier coefficients
fourierofKernel = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform
for fi=1:N
% create sine wave

% sine_wave =2*sin(2*pi* (fi-1) *time + 0);
sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
fourierofKernel(fi) = sum(sine_wave.*data);
end


convFreqDomain=fourierofSignal.*fourierofKernel;

%% Concolution as Filter
lenSin=1000;
timeSin = (0: lenSin-1)/lenSin; 
srate=15;
timeGaus=-2:1/srate: 2;
freqSin=20;
freqGaus=15;
s=6/(2*pi*freqGaus); %s=n/2*pi*f   (The standard deviation of the Gaussian)


sinwave=sin(2*pi* freqSin *timeSin + 0);
gausswave=exp(-(timeGaus.^2)/(2*s^2)); 
figure, plot(gausswave),figure, plot(sinwave)
c=zeros(1, 1000); 
for i=1:939
c(i)=sum(sinwave(i:i+60).*gausswave);
end
figure, plot(sinwave)
hold on, plot(c, 'r-')

%% Fourier transform
%Real EEG data

% sinyal ile ayný zaman uzunluðunda N (sinyalin uzunluðu) tane farklý
% frekans da (1: N freq) sinüs  sinyali yap ve sinyal ile dot product yap.
% Örneðin 5 hz ile dot product sonucu sþnyaldeki 5Hz lik sinyalin genliðini
% verir.
close all; clear all; clc
load('inputData.mat')
data = inputData(1,:); 
N = length(data); 
Fs=1000;
f = Fs/2*linspace(0,1,N/2+1);

% initialize Fourier coefficients
fourier = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform
for fi=1:N
% create sine wave

% sine_wave =2*sin(2*pi* (fi-1) *time + 0);
sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
fourier(fi) = sum(sine_wave.*data);
end

figure, plot(f,2*abs(fourier(1:N/2+1)))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')


% Generated Sum of sine waves
len=1000;
time = (0: len-1)/len; 
data =2*sin(2*pi* 10 *time + 0)+2*sin(2*pi* 30 *time + 0)+2*sin(2*pi* 20 *time + 0);
figure, plot(data);

N = length(data); 
Fs=200;
f = Fs/2*linspace(0,1,N/2+1);

% initialize Fourier coefficients
fourier = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform
for fi=1:N
% create sine wave
 sine_wave =2*sin(2*pi* (fi-1) *time + 0);
% sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
fourier(fi) = sum(sine_wave.*data);
end

figure, plot(abs(fourier (1:50))/N);
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')

%% Inverse Fourier Transform
len=1000;
time = (0: len-1)/len; 
Data = zeros(size(fourier));

for fi=1:len
% create sine wave
 sine_wave =2*sin(2*pi* (fi-1) *time + 0);
% sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
Data(fi,:) = sine_wave*fourier(fi);
end
temporalSignal=mean(Data);
plot(temporalSignal);

 %% Morlet Wavelet
 % sinwave ile gauss wave in point by point çarpýmý
 % s arttýkça gausswave de bell tombullaþýr. Morlet ise sinus a daha çok
 % benzemeye baþlar.
 % a arttýkça gauss un geniþliði artar bu da wavelet in cycle sayýsýný arttýrýr. Bu durum da wavele tin power
 % spectrum da daha keskin frekans bilgisini verir. (freq precision) 
 % This is the parameter that controls the Heisenberg uncertainty principle applied to time-frequency analysis:
 
 % gauss ve sin wave lerin frekanslarýnýn  ayný olmasý gerekir. (peak frequency)
 
 
% srate=100;
%  time=-2:1/srate: 2;
time=linspace(-1,1,100);

sinwave=sin(2*pi* 10 *time + 0);

s=6/(2*pi*10); %s=n/2*pi*f   (The standard deviation of the Gaussian)

gausswave=exp(-(time.^2)/(2*s^2)); 
close all,
plot(gausswave)
figure, plot(sinwave)

morlet=sinwave.*gausswave;
figure, plot(morlet);


%% Morlet Conv as filter
%sinyal ile morlet in convolution ý (conv: sinyal üzerinde kernel i gezdirerek dot product yap)
% wavelet in peak frequency cevresi filtrelenir. (Örn: peak freq: 6Hz ise 4-8 Hz arasý band pass filter gibi)

load ('inputData.mat');
c2=conv(inputData(1,:), morlet);
figure, plot(inputData(1,:)), hold on, plot(c2, 'r');


%% Complex Numbers and Euler's Formula
a=4-8*i;
M=sqrt(4^2+(8)^2)
theta=atan(-8/4)

imageryPart=M*sin(theta);
realPart=M*cos(theta);

M*(cos(theta)+i*sin(theta)) % sin and cos 
M*exp(i*theta) % eular 

%% Complex Morlet
%complex sinus dalgasý (exp(i*2*pi*f*t) ile gauss dalgasýnýn point by point çarpýmý


srate=750;
time=-2:1/srate: 2; time=time(1:3000);
 s=6/(2*pi*10); %s=n/2*pi*f   (The standard deviation of the Gaussian)
A=1/sqrt((s*sqrt(pi))); % frequency band-specific scaling factor.

 %complex sinwave
complexsinwave=exp(i*2*pi*10*time);  % theta=2*pi*f*time (eular' formula)
figure, plot(complexsinwave)

%gauss wave
gausswave=exp(-(time.^2)/(2*s^2)); 
figure, plot(gausswave)

%complex wavelet
complexmorlet=A*(gausswave.*complexsinwave);


%Dot product with data
%fourier gibi 

load('inputData.mat')
data = inputData(1,:); 
N = length(data); 

f =1:Fs/2+1;


% initialize Fourier coefficients
morletDot = zeros(N);
gausswave=exp(-(time.^2)/(2*s^2)); 
% Fourier transform
for fi=1:length(f)

complexsinwave=exp(i*2*pi*(fi-1)*time);
complexmorlet=A*(gausswave.*complexsinwave);
% compute dot product between sine wave and data
morletDot(fi) = sum(complexmorlet.*data);
end

figure, plot(f,abs(morletDot).^2)  % 13.7 D2
xlabel('Time')
ylabel('Frequency')

%% Hilbert Transform
close all
len=1000;
time = (0: len-1)/len; 
cosinwave=cos(2*pi* 10 *time + 0);
figure, plot(cosinwave)
h=hilbert(cosinwave);

himag=imag(h);

hold on, plot(himag, 'r') % 90 derece döndürülmüþ cosinus dalgasý=sinus dalgasý

%% Short-Time Fourier Transform
%Real EEG data
close all; clear all; clc
load('inputData.mat')
data = inputData(1,:); 
data=data(1:250);
H=hann(250);
dataHann=data.*H';
% see removed edge artifact
figure, plot(data), hold on, plot(dataHann, 'r')

N = length(dataHann); 
Fs=250;
f = Fs/2*linspace(0,1,N/2+1);

% initialize Fourier coefficients
fourier = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform
for fi=1:N
% create sine wave

% sine_wave =2*sin(2*pi* (fi-1) *time + 0);
sine_wave = exp(-1i*2*pi*(fi-1).*time);
% compute dot product between sine wave and data
fourier(fi) = sum(sine_wave.*dataHann);
end

figure, plot(f,2*abs(fourier(1:N/2+1)))
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
