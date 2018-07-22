+%% Niruyan Rakulan 214343438
% Introduction: The purpose of this lab was to realize the properties of
% DFT on Matlab, and use filters to remove noise from a sound file.It was dicovered in the lab that
% DFT has hermitian symetry, and noise from music can be removed using filters
% through convolution, and well as multiplication when both signals are in
% frequency domain. The spectogram was used here.

%Materials: PC,Matlab,headphones

%% Q1
%% 1a
clear all;
close all;
x=[2 1 1 0 1];
N=length(x);
M=1000;
X=fft(x,M);
Xmr=fliplr(X);
Xcmr=conj(Xmr);

nlen=length(X);
f=(0:1/nlen:1-1/nlen);

figure
subplot(2,2,1);
stem(f,abs(X));
%*X[M-r]
subplot(2,2,3);
stem(f,abs(Xcmr));

subplot(2,2,2);
stem(f,angle(X));
%X*[M-r]
subplot(2,2,4);
stem(f,angle(Xcmr));
% The frequncy response for both signals are identical. THe DFT was
% conducted with 1000 points to increase the resolution.
%X[r] for the 2-11th point [4.99958549434325 - 0.0439792792922490i,4.99834219009313 - 0.0879404535064596i,4.99627072528201 - 0.131865427913275i,4.99337216287062 - 0.175736128474322i,4.98964799010848 - 0.219534512171746i,4.98510011763912 - 0.263242577318624i,4.97973087835050 - 0.306842373843884i,4.97354302597157 - 0.350316013545391i,4.96653973341556 - 0.393645680304887i,4.95872459087138 - 0.436813640258472i,4.95010160364396 - 0.479802251916369i]
%X*[M-r] for the first 10 points [4.99958549434325 - 0.0439792792922490i,4.99834219009313 - 0.0879404535064596i,4.99627072528201 - 0.131865427913275i,4.99337216287062 - 0.175736128474322i,4.98964799010848 - 0.219534512171746i,4.98510011763912 - 0.263242577318624i,4.97973087835050 - 0.306842373843884i,4.97354302597157 - 0.350316013545391i,4.96653973341556 - 0.393645680304887i,4.95872459087138 - 0.436813640258472i,4.95010160364396 - 0.479802251916369i]
%X[r]=X*[m-r]
%% 1b

figure 
subplot(1,2,1);
stem(f,abs(X));
subplot(1,2,2);
stem(f,abs(Xmr));

% The magnitude response for both signals are identical. THe DFT was
% conducted with 1000 points to increase the resolution.
%% 1c

figure 
subplot(1,2,1);
stem(f,angle(X));
subplot(1,2,2);
stem(f,-angle(Xmr));

% The phase response for both signals are identical. THe DFT was
% conducted with 1000 points to increase the resolution.
%% Q2
clear all;
close all;
x1=[2 1 0 0 0 0 1];
x2=[5 1 1 1 1];
x3=[1 3 3 0 3 3];

N1=length(x1);
N2=length(x2);
N3=length(x3);

f1=(0:1/N1:1-1/N1)*5;
f2=(0:1/N2:1-1/N2)*5;
f3=(0:1/N3:1-1/N3)*5;

X1=fft(x1,N1);
X2=fft(x2,N2);
X3=fft(x3,N3);

figure
subplot(3,2,1);
stem(f1,abs(X1));
subplot(3,2,2);
stem(f1,angle(X1));

subplot(3,2,3);
stem(f2,abs(X2));
subplot(3,2,4);
stem(f2,angle(X2));

subplot(3,2,5);
stem(f3,abs(X3));
subplot(3,2,6);
stem(f3,angle(X3));

%It was realized that the phase response for two of the signals was zero.
%There are not enough samples to give a true representation of the acutal
%phase response.
%% Q2ii

clear all;

x1=[2 1 0 0 0 0 1];
x2=[5 1 1 1 1];
x3=[1 3 3 0 3 3];

N1=16;
N2=16;
N3=16;

f1=(0:1/N1:1-1/N1)*5;
f2=(0:1/N2:1-1/N2)*5;
f3=(0:1/N3:1-1/N3)*5;

X1=fft(x1,N1);
X2=fft(x2,N2);
X3=fft(x3,N3);

figure
subplot(3,2,1);
stem(f1,abs(X1));
subplot(3,2,2);
stem(f1,angle(X1));

subplot(3,2,3);
stem(f2,abs(X2));
subplot(3,2,4);
stem(f2,angle(X2));

subplot(3,2,5);
stem(f3,abs(X3));
subplot(3,2,6);
stem(f3,angle(X3));

% The resolution of the frequency axis was greatly increased; what has
% happened was zeros were padded at the end of the x, increasing the
% resolution of the frequncy response for thier respcetive DFTs.
%% Q3

%% 3.1
clear all;
[y,FS]=audioread('music_noisy.wav');
sound(y,FS);
%The file is of music playing but with 2 distinct tones in the background
%interfering with the music.

%% 3.2
figure;
Y=fft(y);
nlen=length(Y);
f=(0:1/nlen:1-1/nlen)*FS;
stem(f,Y);
%THe noises are the peaks in the DFT.The two noise frequcnies are 1102 2756HZ.
%% 3.3
figure;
load filters;
freqz(h2,1);%notch frequency=0.05pirad/sample*1cycle/2pirad*22050samples/sec=551Hz
figure;
freqz(h4,1);%notch frequency=1102.5Hz
figure;
freqz(h6,1);%notch frequency=1653.75Hz
figure;
freqz(h8,1);%notch frequency=2205Hz
figure;
freqz(h10,1);%notch freqeuncy=2756.25Hz
figure;
freqz(h12,1);%notch freqeuncy=3307.5Hz
figure;
freqz(h14,1);%notch freqeuncy=3858.75Hz
figure;
freqz(h16,1);%notch freqeuncy=4410Hz
%These are bandstop fitlers since one of the frequecncies is attenuated,
%but allows the others.
%% 3.4
%h4 and h10 were chosen becuase it would attuenuate the noise.

%% 3.5
y1=conv(y,h4');
y1=conv(y1,h10');
sound(y1,FS);
% Method 1:Convoltuon of y and h4 and h10.

H4=fft(h4,297702);
H10=fft(h10,297702);
H4=H4';
H10=H10';

Yc=Y.*H4;
Yc=Yc.*H10;
sound(ifft(Yc),FS);
% Method 2: Multiplication of the music signal and the filter signal in
% frequncy domain; converted back to time signal inorder to play music.

%% 3.6
figure;
Y1=fft(y1);
nleny1=length(Y1);
subplot(1,2,1);
f=(0:1/nleny1:1-1/nleny1)*FS;
stem(f,Y1);
%magnitude repsone of method 1; no more noise
nlenyc=length(Yc);
subplot(1,2,2);
f=(0:1/nlenyc:1-1/nlenyc)*FS;
stem(f,Yc);
%magnitude repsone of method 2; no more noise
%% 4.1
clear all;
close all;
[y1,FS1]=audioread('CScale.wav');
[y2,FS2]=audioread('CScaleZ.wav');
figure;
sound(y1,FS1);
sound(y2,FS2);
subplot(1,2,1);
Y1=fft(y1);
nlen1=length(Y1);
f=(0:1/nlen1:1-1/nlen1)*FS1;
stem(f,abs(Y1));

subplot(1,2,2);
Y2=fft(y2);
nlen2=length(Y2);
f=(0:1/nlen2:1-1/nlen2)*FS2;
stem(f,abs(Y2));
% The difference cant be told between the two plots.
%% 4.2
%There are 4002 points in ecah key.
%% 4.3
window_length=4002;
x1=y1(1:4001);
x2=y1(4001:8002);
x3=y1(8002:12003);
x4=y1(12003:16004);
x5=y1(16004:20005);
x6=y1(20005:24006);
x7=y1(24006:28007);
x8=y1(28007:32008);
%magnitude response of key 1
figure;
x1fft=fft(x1);
nlen3=length(x1fft);
f=(0:1/nlen3:1-1/nlen3)*FS1;
stem(f,abs(x1fft));
%magnitude response shows only one peark of 261 Hz

%% 4.4
x1fft=fft(x1,4002);
x2fft=fft(x2,4002);
x3fft=fft(x3,4002);
x4fft=fft(x4,4002);
x5fft=fft(x5,4002);
x6fft=fft(x6,4002);
x7fft=fft(x7,4002);
x8fft=fft(x8,4002);
%% 4.5 Already columns

%% 4.6
spect=[x1fft,x2fft,x3fft,x4fft,x5fft,x6fft,x7fft,x8fft];
%% 4.7
figure
spect_mag=20*log10(abs(spect));
t=(0:window_length:(length(y1)-window_length))/FS1;
f=(1:window_length)*FS1/window_length;
imagesc(t, f, spect_mag);
axis xy
colormap(jet)
colorbar
% It does match the frequncies as time goes on. The frequncies increase
% according to Cscale.
%% Cscalez
window_length=4002;
x1=y2(1:4001);
x2=y2(4001:8002);
x3=y2(8002:12003);
x4=y2(12003:16004);
x5=y2(16004:20005);
x6=y2(20005:24006);
x7=y2(24006:28007);
x8=y2(28007:32008);

x1fft=fft(x1,4002);
x2fft=fft(x2,4002);
x3fft=fft(x3,4002);
x4fft=fft(x4,4002);
x5fft=fft(x5,4002);
x6fft=fft(x6,4002);
x7fft=fft(x7,4002);
x8fft=fft(x8,4002);

spect=[x1fft,x2fft,x3fft,x4fft,x5fft,x6fft,x7fft,x8fft];

figure
spect_mag1=20*log10(abs(spect));
t=(0:window_length:(length(y2)-window_length))/FS2;
f=(1:window_length)*FS2/window_length;
imagesc(t, f, spect_mag1);
axis xy
colormap(jet)
colorbar
%The Cscalez looks correct, backwards of cscale.
%% 5
figure
spectrogram(y1, [], 8000);
figure
spectrogram(y2, [], 8000);
%The x and y axis are switched. The frequncy axis is normailized.

%% Conclusion
% The lab went as suspected. THe lab taught use how to normalize
% frequncies, the Hermitian symetry of DFT, how changing the number of
% points of DFT enchances the resolution of DFT, and the purpose of the
% spectogram.