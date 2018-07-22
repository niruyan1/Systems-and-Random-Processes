%% Niruyan Rakulan 214343438
%Introduction: The purpose of the project was to understand FIR,IFIR,
%and FRM filters. Our knowledge of filters was used to clean up a ECG
%signal, which is vital for doctors to diagnose heart conditions.
%Equipment: MATLAB, PC
%% Part 1
clear all;
close all;
samp_freq = 8000;% sample frequency
pb= 3600;% passband of overall filter
sb = 3200;% stopband of overall filter
sr = 0.0001;% stopband ripple
pr = 0.01; % Passband ripple
npb = pb/samp_freq;%normalized passband of overall filter
nsb = sb/samp_freq;%normalized stopband of overall filter
B=npb-nsb;%difference between passband and stopband of normalized filter

%To ensure M(interpolation factor) is even; wont work if odd
%a)M=2

if(mod(floor(1/(1-(npb+nsb)+sqrt(B))),2)==0)
    M=floor(1/(1-(npb+nsb)+sqrt(B)));%M is initially even, stays even
else
    M=ceil(1/(1-(npb+nsb)+sqrt(B)));%M is initially odd, becomes even
end

%b)Hapb=0.1,Hasb=0.2,Hmaskingpb=0.45,Hmaskingsb=0.1
if((round(M*npb)+0.5-M*npb)<=0.5)
    Hapb=0.5-(round(M*npb)+0.5-M*npb);
    Hasb=0.5-(round(M*nsb)+0.5-M*nsb);
    Hmpb=npb;
    Hmsb=(round(M*npb)+0.5-1+(0.5-Hasb))/M;
    
else
    Hapb=round(M*npb)-M*npb;
    Hasb=round(M*nsb)-M*nsb;
    Hmpb=npb;
    Hmsb= Hasb/M;
    
end

%c)La=34,Lm=11
%length of shaping filter
La=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(Hapb-Hasb))+1);%length of bandedge shaping filter
%length of masking filter
Lm=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(Hmpb-Hmsb))+1);%length of masking filter

%d)Overall filter using Ifir.Length(Ha^m)=73,Length(HMa)=11; same value.
%you have to take into consideration that length of Ha^m is almost twice that of
%Ha

%use 'ifir'
[ham,hma,d]=ifir(M,'high',[nsb*2 npb*2],[sr,pr]);
%top branch
branch1=conv(ham,hma);
%frequency repsonse of ham
[Ham,w]=freqz(ham,1);
%frequency response of hma
[Hma,w]=freqz(hma,1);
[Branch1,w]=freqz(branch1,1);
[Branch2,w]=freqz(d,1);
%introduce delay to get high pass
Hovr=Branch1+Branch2;

%plot Ha^M frequncy response
figure;
%magnitued response
plot((samp_freq).*w./(2*pi),abs(Ham));
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('Ha(z^M)(Magnitude Response)Linear');
grid on;
%plot Ha^M using freqz
figure;
freqz(ham,1);
title('Ha(z^M)(Frequency Response) using Freqz');
grid on;

%plot HMa(masking filter) freq. response
figure;
%magnitude response
plot((samp_freq).*w./(2*pi),abs(Hma));
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('HMa(z)(Magnitude Response)Linear');
grid on;
%plot Ha^M using freqz
figure;
freqz(hma,1);
title('HMa(z)(Frequency Response)Using Freqz');
grid on;

%plot frequency response for overall filter
figure;
%magnitude response
hold on;
plot((samp_freq).*w./(2*pi),abs(Hovr));
%make sure parameters are met
plot((samp_freq).*w./(2*pi),ones(size((samp_freq).*w./(2*pi))) * (1+pr));
plot((samp_freq).*w./(2*pi),ones(size((samp_freq).*w./(2*pi))) * (1-pr));
plot((samp_freq).*w./(2*pi),ones(size((samp_freq).*w./(2*pi))) * sr);
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('Hoverall (Magnitude Response)Linear');
grid on;
hold off;

%show the passband meets parameters
figure
ifirmax=imread('q1maxofifirfilter.png');
imshow(ifirmax);
title('Pass band ripple less than 0.01(IFIR Matlab function)');
%show stopband meets parameters
figure
ifirmin=imread('q1minofifirfilter.png');
imshow(ifirmin);
title('Stop band ripple less than 0.0001(IFIR Matlab function)');
%satisfies parameters
lenh=length(ham);%length of bandedge using 'ifir'
leng=length(hma);%length of masking filter using 'ifir'


%a)M=2
%b)Shaping filter passband(Hapb)=0.1,Shaping filter stop band(Hasb)=0.2,
%Masking filter passband (Hmpb)=0.45,Masking filter stopband(Hmsb)=0.1
%c)Length of Ha=34,Masking=11
%d)Length of interpolated=73,Length of Masking=11; is pretty much the same.
%Length of interpolated is 2 times length of Ha since interpoalted by 2
%e)(shown above)

%% 2
clear all;
close all;
%overall filter passband
pb=0.194;
%overall filter stopband
sb=0.2;
%passband ripple
pr=0.01;
%stopband ripple
sr=0.001;
%difference between stopband and passband
B=sb-pb;
%a)Mopt=6
M=round(1/(2*sqrt(sb-pb)));
samp_freq=1;

%initialize Mmin and Lmin for M=6
m=floor(pb*M);
theta=pb*M-m;
phi=sb*M-m;
fpMa=pb;
fsMa=(m+1-phi)/M;
fpMc=(m-theta)/M;
fsMc=sb;
La=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(theta-phi))+1);
Lma=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMa-fsMa))+1);
Lmc=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMc-fsMc))+1);
Lt=La+Lma+Lmc;
Lmin=Lt;
Min=M;
casezmin='A';

%find Lmin for M-2 to M+2
for Mopt=M-2:1:M+2
    %for case a
    % m value for case A
    ma=floor(pb*Mopt);
    %pass and stop band for bandedge filter for case A
    thetaa=pb*Mopt-ma;
    phia=sb*Mopt-ma;
    %masking filter passband for case A
    fpMaa=pb;
    %masking filter stopband for case A
    fsMaa=(ma+1-phia)/Mopt;
    %masking compliment filter passband for case A
    fpMca=(ma-thetaa)/Mopt;
    %masking compliment filter stopband for case A
    fsMca=sb;
    %lenghts for each filter in case A
    Laa=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(thetaa-phia))+1);
    Lmaa=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMaa-fsMaa))+1);
    Lmca=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMca-fsMca))+1);
    Lta=Laa+Lmaa+Lmca;
    caseza='A';
    
    %for case b
    mb=ceil(sb*Mopt);
    %passband for bandedge shaping filter
    thetab=mb-sb*Mopt;
    phib=mb-pb*Mopt;
    %passband for masking filter case B
    fpMab=(mb-1+phib)/Mopt;
    %stopband for masking filter case B
    fsMab=sb;
    %passband for compliment masking filter
    fpMcb=pb;
    %stopband for compliment masking filter
    fsMcb=(mb+thetab)/Mopt;
    %lenghths of each filter
    Lab=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(thetab-phib))+1);
    Lmab=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMab-fsMab))+1);
    Lmcb=ceil(( -20* log10(sqrt(sr*pr))-13)/(14.6*abs(fpMcb-fsMcb))+1);
    Ltb=Lab+Lmab+Lmcb;
    
    casezb='B';
    
    %if A is not suitable
    if((thetaa>=0.5||phia>=0.5)||(thetaa<=0||phia<=0))
        %if A is not suitable, then case assigned 'x' and all edges=0
        ma=0;
        thetaa=0;
        phia=0;
        fpMaa=0;
        fsMaa=0;
        fpMca=0;
        fsMca=0;
        Laa=0;
        Lmaa=0;
        Lmca=0;
        Lta=0;
        caseza='x';
        
    end
    %if B is not suitable
    if((thetab>=0.5||phib>=0.5)||(thetab<=0||phib<=0))
        %if B is not suitable, then case assigned 'x' and all edges=0
        mb=0;
        thetab=0;
        phib=0;
        fpMab=0;
        fsMab=0;
        fpMcb=0;
        fsMcb=0;
        Lab=0;
        Lmab=0;
        Lmcb=0;
        Ltb=0;
        
        casezb='x';
    end
    %print out values for each M value from M-2 to M+2 for both cases
    fprintf('%c %f %f %f %f %f %f %f %f %f %f %f\n', caseza,Mopt,thetaa,phia,fpMaa,fsMaa,fpMca,fsMca,Laa,Lmaa,Lmca,Lta);
    fprintf('%c %f %f %f %f %f %f %f %f %f %f %f\n', casezb,Mopt,thetab,phib,fpMab,fsMab,fpMcb,fsMcb,Lab,Lmab,Lmcb,Ltb);
    %if Lmin found, replace current Lmin
    if(Lta<Lmin&&Lta>0)
        Lmin=Lta;
        Min=Mopt;
        casezmin='A';
    end
    if(Ltb<Lmin&&Ltb>0)
        Lmin=Ltb;
        Min=Mopt;
        casezmin='B';
    end
end
%{
B 4.000000 0.200000 0.224000 0.056000 0.200000 0.194000 0.300000 107.000000 19.000000 25.000000 151.000000
A 6.000000 0.164000 0.200000 0.194000 0.300000 0.139333 0.200000 72.000000 25.000000 43.000000 140.000000
A 7.000000 0.358000 0.400000 0.194000 0.228571 0.091714 0.200000 62.000000 75.000000 25.000000 162.000000
B 8.000000 0.400000 0.448000 0.181000 0.200000 0.194000 0.300000 54.000000 135.000000 25.000000 214.000000
everything else is invalid
%}

%using the new found Mopt, find bandedges
%cMopt=Mmin=6
%using the newly found Moptmin, find the bandedges, if case A
if(casezmin=='A')
    %m for final bandedge shaping filter Case A
    m=floor(pb*Min);
    %passband for final bandedge shaping filter Case A
    theta=pb*Min-m;
    %stopband for final bandedge shaping filter Case A
    phia=sb*Min-m;
    %passband for final masking filter to Case A
    fpMa=pb;
    %stopband for final masking filter to Case A
    fsMa=(m+1-phi)/Min;
    %passband for complement masking filter Case A
    fpMc=(m-theta)/Min;
    %stopband for complement masking filter Case A
    fsMc=sb;
end

%using the newly found Moptmin, find the bandedges, if case B
if(casezmin=='B')
    %m for final bandedge shaping filter Case B
    m=ceil(sb*Min);
    %passband for final bandedge shaping filter Case B
    theta=m-sb*Min;
    phi=m-pb*Min;
    %passband for masking filter case B
    fpMa=(m-1+phi)/Min;
    %stopband for masking filter case B
    fsMa=sb;
    %passband for compliment masking filter case B
    fpMc=pb;
    %stopband for compliment masking fitler case B
    fsMc=(m+theta)/Min;
end

%d and %e
%make a filter for bandedge shaping filter
[Ni1,Fi1,Ai1,W1] = firpmord([theta*samp_freq, phi*samp_freq],[1,0],0.85.*[pr,sr],samp_freq);
ha=firpm(Ni1+9,Fi1,Ai1,W1);
%interpolate ha by M
ham=upsample(ha,M);
ham=ham(1:493);
%plot frequency response of Ha, both linear and in dB
[Ha,wori]=freqz(ha,1);
figure;
hold on;
plot(samp_freq.*wori./(2*pi),abs(Ha));
plot((samp_freq).*wori./(2*pi),ones(size((samp_freq).*wori./(2*pi))) * (1+0.85*pr));
plot((samp_freq).*wori./(2*pi),ones(size((samp_freq).*wori./(2*pi))) * (1-0.85*pr));
plot((samp_freq).*wori./(2*pi),ones(size((samp_freq).*wori./(2*pi))) * 0.85*sr);
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('Ha(z)(Magnitude Response)Linear');
hold off;
grid on;
figure
freqz(ha,1);
title('Ha(z)(Frequency Response)Freqz');


%plot frequency response of interpolated filter
[Ham,w1]=freqz(ham,1);
Hamlen=length(Ham);
%plot in frequency domain Ha^m
figure;
hold on;
plot(samp_freq.*w1./(2*pi),abs(Ham));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1+0.85*pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1-0.85*pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * 0.85*sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Ha(z^M)(Magnitude Response)');
grid on;
hold off;
%using freqz
figure;
freqz(ham,1);
title('Ha(z^M)(Frequency Response)Freqz');

%make a filter for Masking filter
[Ni2,Fi2,Ai2,W2] = firpmord([fpMa*samp_freq, fsMa*samp_freq],[1,0],0.85.*[pr,sr],samp_freq);
hma=firpm(Ni2+9,Fi2,Ai2,W2);
%pad with zeros
hma=[zeros(1,9),hma,zeros(1,9)];
%maskign filter in frequency domain
[Hma,w2]=freqz(hma,1);
%plot frequency response of masking fitler
figure;
hold on;
plot(samp_freq.*w2./(2*pi),(abs(Hma)));
plot((samp_freq).*w2./(2*pi),ones(size((samp_freq).*w2./(2*pi))) * (1+0.85*pr));
plot((samp_freq).*w2./(2*pi),ones(size((samp_freq).*w2./(2*pi))) * (1-0.85*pr));
plot((samp_freq).*w2./(2*pi),ones(size((samp_freq).*w2./(2*pi))) * 0.85*sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('HMa(z)(Magnitude Response)');
grid on;
hold off;
%plot phase response
figure;
freqz(hma,1);
title('HMa(z)(Magnitude Response)Freqz');

%plot compliment masking filter
[Ni3,Fi3,Ai3,W3] = firpmord([fpMc*samp_freq, fsMc*samp_freq],[1,0],0.85.*[pr,sr],samp_freq);
hmc=firpm(Ni3+8,Fi3,Ai3,W3);
%compliment masking filter in frequency domain
[Hmc,w3]=freqz(hmc,1);
figure;
%mangitude response
hold on;
plot(samp_freq.*w3./(2*pi),(abs(Hmc)));
plot((samp_freq).*w3./(2*pi),ones(size((samp_freq).*w3./(2*pi))) * (1+0.85*pr));
plot((samp_freq).*w3./(2*pi),ones(size((samp_freq).*w3./(2*pi))) * (1-0.85*pr));
plot((samp_freq).*w3./(2*pi),ones(size((samp_freq).*w3./(2*pi))) * 0.85*sr);
hold off;
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('HMc(z)(Magnitude Response)Linear');
grid on;
figure;
freqz(hmc,1);
title('HMc(z)(Magnitude Response)Freqz');

%plot compliment filter
Hcm=exp(-1i.*w1.*(length(ha)-1).*Min./2)-Ham;
figure
hold on;
plot(samp_freq.*w1./(2*pi),(abs(Hcm)));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1+0.85*pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1-0.85*pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * 0.85*sr);
%because it was compliment, the pass ripple will be on the stop band for
%compliment and vice versa
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Hc(z^M)(Magnitude Response) Linear');
grid on;
hold off;
axis([0 0.5 0 1.2]);

%overall filter
%lengths of filter changed to meet specifications
%Ha=83,HMa=52,HMc=52

%f and %g
Hovr=Ham.*Hma+Hcm.*Hmc;
figure
hold on;
plot(samp_freq.*w1./(2*pi),(abs(Hovr)));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1+pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * (1-pr));
plot((samp_freq).*w1./(2*pi),ones(size((samp_freq).*w1./(2*pi))) * sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude Response')
title('Hovr(Magnitude Response');
grid on;
hold off;
figure
ifirmin=imread('q2ovrhigh.png');
imshow(ifirmin);
title('Pass band ripple less than 0.01');
figure
ifirmin=imread('q2ovrlow.png');
imshow(ifirmin);
title('Stop band ripple less than 0.001');

%a)6
%b)
%{
B 4.000000 0.200000 0.224000 0.056000 0.200000 0.194000 0.300000 107.000000 19.000000 25.000000 151.000000
A 6.000000 0.164000 0.200000 0.194000 0.300000 0.139333 0.200000 72.000000 25.000000 43.000000 140.000000
A 7.000000 0.358000 0.400000 0.194000 0.228571 0.091714 0.200000 62.000000 75.000000 25.000000 162.000000
B 8.000000 0.400000 0.448000 0.181000 0.200000 0.194000 0.300000 54.000000 135.000000 25.000000 214.000000
everything else is invalid
%}
%c)6 Case A considered optimal interpolation factor, length of 140
%d)shown above
%e)shown above
%f)shown above
%g)shown above; filter lengths changed:%Ha=83,HMa=52,HMc=52
%% part3
clear all;
close all;
samp_freq = 360;% samplefreq
pr = 0.01; % Passband ripple
sr = 0.001;% stopband ripple
%filter1(highpass fitler)
%passband for highpass filter 1
sb1=0.17;
%stopband for filter 1
pb1=1;
%normalized stopband and passband filter
nsb1=sb1/samp_freq;
npb1=pb1/samp_freq;
F1=[sb1,pb1];
A1=[0,1];

%use of firpm function for highpass(0.17 and 1)
[Ni1,Fi1,Ai1,Wi1]=firpmord(F1,A1,[sr,pr],samp_freq);
h1=firpm(Ni1+8,Fi1,Ai1,Wi1);
[Hh,wh]=freqz(h1,1);
%plot freq. response of filter 1
figure
hold on;
plot(samp_freq.*wh./(2*pi),(abs(Hh)));
plot((samp_freq).*wh./(2*pi),ones(size((samp_freq).*wh./(2*pi))) * (1+pr));
plot((samp_freq).*wh./(2*pi),ones(size((samp_freq).*wh./(2*pi))) * (1-pr));
plot((samp_freq).*wh./(2*pi),ones(size((samp_freq).*wh./(2*pi))) * sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Highpass Filter(Magnitude Response) Linear');
hold off;
figure
freqz(h1,1);
title('Filter 1, Highpass Filter')



%filter2(low pass filter)
%stop band for filter 2
sb2=40;
%passband for filter 2
pb2=54;
F2=[sb2,pb2];
A2=[1,0];
[Ni2,Fi2,Ai2,Wi2]=firpmord(F2,A2,[pr,sr],samp_freq);
h2=firpm(Ni2+20,Fi2,Ai2,Wi2);
[Hl,wl]=freqz(h2,1);
figure
hold on;
plot(samp_freq.*wl./(2*pi),(abs(Hl)));
plot((samp_freq).*wl./(2*pi),ones(size((samp_freq).*wl./(2*pi))) * (1+pr));
plot((samp_freq).*wl./(2*pi),ones(size((samp_freq).*wl./(2*pi))) * (1-pr));
plot((samp_freq).*wl./(2*pi),ones(size((samp_freq).*wl./(2*pi))) * sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Lowpass Filter(Magnitude Response) Linear');
hold off;
figure
freqz(h2,1);
title('Filter 2, Lowpass Filter');




load('ECG_signal.mat');

%Overall filter
hovr=conv(h1,h2);
[Hovr,wovr]=freqz(hovr,1);
figure
hold on;
plot(samp_freq.*wovr./(2*pi),(abs(Hovr)));
plot((samp_freq).*wovr./(2*pi),ones(size((samp_freq).*wovr./(2*pi))) * (1+pr));
plot((samp_freq).*wovr./(2*pi),ones(size((samp_freq).*wovr./(2*pi))) * (1-pr));
plot((samp_freq).*wovr./(2*pi),ones(size((samp_freq).*wovr./(2*pi))) * sr);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Overall Filter(Magnitude Response) Linear');
hold off;
%frequency response of overall filter
figure
freqz(hovr,1);
title('Overall Filter')


figure
ifirmax=imread('q3passband.png');
imshow(ifirmax);
title('Pass band ripple less than 0.01');
%show stopband meets parameters
figure
ifirmin=imread('q3stopband.png');
imshow(ifirmin);
title('Stop band ripple less than 0.001');
%satisfies parameters

%plot freq.reposne of original signal
filtecg=filter(hovr,1,ecg');
figure
freqz(ecg,1);
title('Freq. Response of Original Signal');
%plot freq.reponse of filtered signal
figure
freqz(filtecg,1);
title('Freq. Response of Filtered Signal');

%plot oringal singal in time domain
figure
subplot(1,2,1);
plot(t,ecg');
xlabel('Time')
ylabel('Magnitude')
title('Original Ecg');
grid on;
%plot filtered signal in time domain
subplot(1,2,2);
plot(t,filtecg);
xlabel('Time')
ylabel('Magnitude')
title('Filtered Ecg');
grid on;

%length of h1 filter big, probably better to use complementary ifir low
%pass filter

%% %% Conclusion
%The project went mostly as planned. Matlab uses a sample frequency of 2, and
%modifications had to be made inorder to satisfy requirements. Lengths of
%the filters might have been increased inorder to meet ripple parameters
%The project demonstrated the uses of a filter