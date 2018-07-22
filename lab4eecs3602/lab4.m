%% Niruyan Rakulan 214343438
%Introduction: The purpose of this lab was to create a filter using Matlabs
%firpm function. The filter was able to filter the high end frequencies in
%the sound file, and only allowed the low frequencies.

%Materials: PC, Matlab, Headphones
%%  Q1
%% 1)
clear all;
samp_freq = 22050;
pB = 3000;% passband
sB = 5000;% stopband
sr = 0.001;% stopband ripple
pr = 0.01; % Passband ripple
fp = pB/samp_freq;%normalized pB
fs = sB/samp_freq;%normalized sB
%Use of Equation 1
Lp =( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fs-fp))+1; %length
%Lp=28.94

F = [pB, sB];
A = [1,0];
Dev = [pr, sr];
%Use of MATLAB function

L = firpmord(F,A,Dev,samp_freq );
%The two methods did not yield the same answer; Lp(Equation
%1)=28.94,L(Matlab)=27

%% 2)
%FIR filter using equation 1
figure
[q,Fi,Ai,W]=firpmord(F,A,Dev,samp_freq); 
h1=firpm(ceil(Lp),Fi,Ai,W); %filter coeffcient
freqz(h1,1);

figure 
lpbwq1=imread('lowpassbandequation1.png');
imshow(lpbwq1);
%Low pass band for Equation 1. Does meet requirement since lower than
%60dB(-20log(0.001)=60). 

figure
hpbwq1=imread('highpassbandequation1.png');
imshow(hpbwq1);
%High pass band for Equation 1. Does meet requirement since lower than
%0.0864dB(-20log(0.01)=0.0864). 
% The FIR filter design does meet the specifications because the ripples
% are lower than needed.
%% 3)
%FIR filter using matlab equation
figure
[N,Fi,Ai,W]=firpmord(F,A,Dev,samp_freq); 
h2=firpm(N,Fi,Ai,W); %filter coeffcient
freqz(h2,1); %plot

figure 
lpbwm=imread('lowpassbandmatlab.png');
imshow(lpbwm);
%Low pass band for Matlab function. Doesnt meet requirement since higher than
%60dB(-20log(0.001)=60). 

figure
hpbwm=imread('highpassbandmatlab.png');
imshow(hpbwm);
%High pass band for Matlab function. Doesnt meet requirement since higher than
%0.0864dB(-20log(0.01)=0.0864). 
%This filter does not meet the specifications because the ripples are
%higher than intended.
%% 4

%Equation 1 does meet the requirement(assume use ceiling of Lm) since ripples less than intended. Matlab funciton does not meet requirements since ripples greater than intended. To improve, increase the length.
%% 5
figure;
[N,Fi,Ai,W]=firpmord(F,A,Dev,samp_freq); 
h2=firpm(N+2,Fi,Ai,W); %filter coeffcient
freqz(h2,1);
%Length of Matlab function increased by 2.
figure;
lpbwq1=imread('lowpassbandequation1.png');
imshow(lpbwq1);
%Low pass band for New signal. Does meet requirement since lower than
%60dB(-20log(0.001)=60). 
figure
hpbwq1=imread('highpassbandequation1.png');
imshow(hpbwq1);
%High pass band for New Signal. Does meet requirement since lower than
%0.0864dB(-20log(0.01)=0.0864). 
%% 6
figure;
[y,FS] = audioread('music.wav');
sound(y,FS);
Y = fft(y);
m = abs(Y);
n = length(Y);
f = (0:1/n:1-1/n)*FS;
plot(f,m);

%6aconvolition
figure;
yhcon = conv(y,h2);
sound(yhcon,FS);
plot(abs(fft(yhcon)));
% Got rid of the high frequncies, and kept the low frequencies intact.

%bfft
figure;
H2=fft(h2,297702);
YH2Mul= Y.*H2';
sound(ifft(YH2Mul),FS);
plot(abs(YH2Mul));
%Got rid of the high frequncies, and kept the low frequencies intact.

%c
figure;
m1=filter(h2,1,y);
sound(m1,FS);
plot(abs(fft(m1)));
%Got rid of the high frequncies, and kept the low frequencies intact.

%% Conclusion
%The lab went as expected. The lab taught us how to make filters and apply
%them in matlab.