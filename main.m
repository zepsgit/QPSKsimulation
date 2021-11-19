clear;

%%%%%%%%%%%%%%%%%%%%%

%%% Initialize parameters

%%%%%%%%%%%%%%%%%%%%

T=1;                          % the width of baseband, namely the frequecy of baseband

fc=10/T;                      % carrier frequency

ml=2;                         % for separating odd and even positioning number

nb=50;                       % the number of bit to transmit

delta_T=T/200;                % sampling period

fs=1/delta_T;                 % sampling frequency

SNR=0.07;                        % signal noise ratio

t=0:delta_T:nb*T-delta_T;     % the time slot 

N=length(t);                  % the number of sampling points 
% root raised cosine parameters
rf=0.35;%roll off factor
span=1;
sps=fs;
rcos = rcosdesign(rf,span,sps,'sqrt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caculate the energy of pulse shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eg=sum(rcos(1:fs).*rcos(1:fs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%   modulation part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% root raised cosine pulse shape
figure;
plot(rcos);
title("root raised cosine pulse shape");

% generate baseband signal

data=randn(1,nb)>0.5;  % generate 0 or 1 as our data
data_=zeros(1,nb/delta_T);
for q=1:nb

    data_((q-1)/delta_T+1:q/delta_T)=data(q); % transfer binary code to waveform

end
figure;
plot(data_);
ylim([-0.5, 1.5]);
title("original data waveform")
datanrz=data.*2-1;            % transfer to polar code

data1=zeros(1,nb/delta_T);    % generate a matrix for storing data

for q=1:nb

    data1((q-1)/delta_T+1:q/delta_T)=datanrz(q)*rcos(1:fs); % transfer polar code to waveform

end

%plot polar code
%{
figure
plot(data1);
%ylim([-1.5,1.5]);
title("trasmitted polar code waveform");
%}
% baseband data to waveform
data0=zeros(1,nb/delta_T);   

for q=1:nb

    data0((q-1)/delta_T+1:q/delta_T)=data(q)*rcos(1:fs); 

end

%plot baseband data waveform
%{
figure;clear;

%%%%%%%%%%%%%%%%%%%%%

%%% Initialize parameters

%%%%%%%%%%%%%%%%%%%%

T=1;                          % the width of baseband, namely the frequecy of baseband

fc=10/T;                      % carrier frequency

ml=2;                         % for separating odd and even positioning number

nb=10;                       % the number of bit to transmit

delta_T=T/200;                % sampling period

fs=1/delta_T;                 % sampling frequency

SNR=0.07;                        % signal noise ratio

t=0:delta_T:nb*T-delta_T;     % the time slot 

N=length(t);                  % the number of sampling points 
% root raised cosine parameters
rf=0.35;%roll off factor
span=1;
sps=fs;
rcos = rcosdesign(rf,span,sps,'sqrt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caculate the energy of pulse shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eg=sum(rcos(1:fs).*rcos(1:fs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%   modulation part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% root raised cosine pulse shape
figure;
plot(rcos);
title("root raised cosine pulse shape");

% generate baseband signal

data=randn(1,nb)>0.5;  % generate 0 or 1 as our data
data_=zeros(1,nb/delta_T);
for q=1:nb

    data_((q-1)/delta_T+1:q/delta_T)=data(q); % transfer binary code to waveform

end
figure;
plot(data_);
ylim([-0.5, 1.5]);
title("original data waveform")
datanrz=data.*2-1;            % transfer to polar code

data1=zeros(1,nb/delta_T);    % generate a matrix for storing data

for q=1:nb

    data1((q-1)/delta_T+1:q/delta_T)=datanrz(q)*rcos(1:fs); % transfer polar code to waveform

end

%plot polar code
%{
figure
plot(data1);
%ylim([-1.5,1.5]);
title("trasmitted polar code waveform");
%}
% baseband data to waveform
data0=zeros(1,nb/delta_T);   

for q=1:nb

    data0((q-1)/delta_T+1:q/delta_T)=data(q)*rcos(1:fs); 

end

%plot baseband data waveform
%{
figure;
plot(data0);
%ylim([-1.5,1.5]);
title("baseband data waveform");
% transmitted signal
%}
data2=abs(fft(data1));
figure;
plot(data2);
title("baseband spectrum");
plot(data0);
%ylim([-1.5,1.5]);
title("baseband data waveform");
% transmitted signal
%}
data2=abs(fft(data1));
figure;
plot(data2);
title("baseband spectrum");

% serial to parrallel
idata=datanrz(1:ml:(nb-1));   % separate odd and even so m1=2

qdata=datanrz(2:ml:nb);

% QPSK modulation

ich=zeros(1,nb/delta_T/2);    % generate a matrix of 1*nb/delta_T/2 for storing odd and even position data

for i=1:nb/2

    ich((i-1)/delta_T+1:i/delta_T)=idata(i)*rcos(1:fs);

end

for ii=1:N/2

    a(ii)=cos(2*pi*fc*t(ii));

end

idata1=ich.*a;                % obtain i data = Am*gt*cos
%plot i data waveform
%{
figure;
plot(idata1);
title("i data waveform");
%}
qch=zeros(1,nb/2/delta_T);

for j1=1:nb/2

    qch((j1-1)/delta_T+1:j1/delta_T)=qdata(j1)*rcos(1:fs);

end

for jj=1:N/2

    b(jj)=sin(2*pi*fc*t(jj));

end

qdata1=qch.*b;               % modulated q data
%plot q data
%{
figure;
plot(qdata1);
title("q data waveform");
%}

s=idata1+qdata1;             % combine i and q
%plot modulated signal
figure;
plot(s );
title("modulated waveform");
%plot modulated signal spectrum
ss=abs(fft(s));             
figure;
plot(ss);
title("i plus q signal spectrum");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rayleigh channel

ray_ich=raylrnd(0.1,1,nb/2/delta_T);
%{
figure;
plot(ray_ich);
title("Rayleigh noise waveform");
%}
ray_qch=raylrnd(0.1,1,nb/2/delta_T);

Ray_idata=idata1+ray_ich;

Ray_qdata=qdata1+ray_qch;

Ray_s=Ray_idata+Ray_qdata;
figure;
plot(Ray_s);
title("transmitted signal after Rayleigh channel");
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian channel

s1=awgn(s,SNR,'measured','linear');              % signal after Gaussian channel
figure;
plot(s1);
title("transmitted signal after Gaussian channel");
s11=abs(fft(s1));  
%{
figure;
plot(s11);
title("transmitted signal spectrum after Gaussian");

%}
s111=s1-s;                   % the waveform of Gaussian noise
%{
figure;
plot(s111);
title("Gaussian noise waveform");
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% filter the received signal
fl=8;%low frequency of filter
fh=15;%hign frequency of filter
% fs--sampling frequency
% will plot the filtered signal automatically
bandpass(s1,[fl fh],fs);
filtered_s=bandpass(s1,[fl fh],fs);
%{
figure;
plot(abs(fft(filtered_s)));
title("filtered transmitted signal spectrum");
%}
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demodulate part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idata2=zeros(1,nb/2/delta_T);
idata2_=(2/Eg)*s1.*a;               % a is cos
for j1=1:nb/2

    idata2((j1-1)/delta_T+1:j1/delta_T)=idata2_((j1-1)/delta_T+1:j1/delta_T).*rcos(1:fs);

end
qdata2=zeros(1,nb/2/delta_T);
qdata2_=(2/Eg)*s1.*b;                % b is sin
for j1=1:nb/2

    qdata2((j1-1)/delta_T+1:j1/delta_T)=qdata2_((j1-1)/delta_T+1:j1/delta_T).*rcos(1:fs);

end

figure;
plot(idata2);
title("demodulated i data waveform");
%{
figure;
plot(qdata2);
title("demodulated q data waveform");
%} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idata4=zeros(1,nb/2);

qdata4=zeros(1,nb/2);

for n=1:nb/2

    Awgn_ichsum(n)=sum(idata2((n-1)/delta_T+1:n/delta_T));

    if Awgn_ichsum(n)>=0

        idata4(n)=1;

    else idata4(n)=0;

    end

    Awgn_qchsum(n)=sum(qdata2((n-1)/delta_T+1:n/delta_T));

    if Awgn_qchsum(n)>=0

        qdata4(n)=1;

    else qdata4(n)=0;

    end

end
figure;
plot(Awgn_ichsum,Awgn_qchsum,'rx')
title("received signal constellation");
% store determined data
%{
figure;
plot(idata4*2-1,qdata4*2-1,'ro');
axis([-2 2 -2 2]);
title("constellation after decision");
%}
demodata=zeros(1,nb);

demodata(1:ml:(nb-1))=idata4; % for odd positioning number

demodata(2:ml:nb)=qdata4;     % for even positioning number

%transform waveform with 1 width for display

demodata1=zeros(1,nb/delta_T);   

for q=1:nb

    demodata1((q-1)/delta_T+1:q/delta_T)=demodata(q); 

end
figure;
plot(demodata1);
ylim([-0.5 1.5])
title("demodulated received signal");

% the number of error

% abs(demodata-data) if the difference is 0 then no error, error happend
% when it is 1

Awgn_num_BER=sum(abs(demodata-data));

%%%%%%%%%%%%%%%%%%%

% Rayleigh demodulation

Ray_idata2_=(2/Eg)*Ray_s.*a;                

Ray_qdata2_=(2/Eg)*Ray_s.*b;

Ray_idata2=zeros(1,nb/2/delta_T);
Ray_qdata2=zeros(1,nb/2/delta_T);
for j1=1:nb/2

    Ray_idata2((j1-1)/delta_T+1:j1/delta_T)=Ray_idata2_((j1-1)/delta_T+1:j1/delta_T).*rcos(1:fs);
    Ray_qdata2((j1-1)/delta_T+1:j1/delta_T)=Ray_qdata2_((j1-1)/delta_T+1:j1/delta_T).*rcos(1:fs);
end

% decision

Ray_idata4=zeros(1,nb/2);

Ray_qdata4=zeros(1,nb/2);

for n=1:nb/2

    Ray_ichsum(n)=sum(Ray_idata2((n-1)/delta_T+1:n/delta_T));

    if Ray_ichsum(n)>=0

        Ray_idata4(n)=1;

    else Ray_idata4(n)=0;

    end

    Ray_qchsum(n)=sum(Ray_qdata2((n-1)/delta_T+1:n/delta_T));

    if Ray_qchsum(n)>=0

        Ray_qdata4(n)=1;

    else Ray_qdata4(n)=0;

    end

end

 

% storing the determined number

Ray_demodata=zeros(1,nb);

Ray_demodata(1:ml:(nb-1))=Ray_idata4; 

Ray_demodata(2:ml:nb)=Ray_qdata4;    


%tranform to waveform for display
Ray_demodata1=zeros(1,nb/delta_T);    

for q=1:nb

    Ray_demodata1((q-1)/delta_T+1:q/delta_T)=Ray_demodata(q); %transform to waveform

end

% accumulated error bit

Ray_num_BER=sum(abs(Ray_demodata-data));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%                BER caculator

%%     involke cm_sm32() and cm_sm33(）

%% cm_sm32() for Rayleigh

%%      cm_sm33() for Gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNRindB1=0:1:6;

SNRindB2=0:0.1:6;

% Rayleigh fading channel

for i=1:length(SNRindB1)

    [pb,ps]=cm_sm32(SNRindB1(i));        

    smld_bit_ray_err_prb(i)=pb;

    smld_symbol_ray_err_prb(i)=ps;

    disp([ps,pb]);

end;

% Gaussian

for i=1:length(SNRindB1),

    [pb1,ps1]=cm_sm33(SNRindB1(i));            

    smld_bit_awgn_err_prb(i)=pb1;

    smld_symbol_awgn_err_prb(i)=ps1;

    disp([ps1,pb1]);

end;

% theoretical curve

for i=1:length(SNRindB2),

    SNR=exp(SNRindB2(i)*log(10)/10);                 % linear SNR

    theo_err_awgn_prb(i)=0.5*erfc(sqrt(SNR));        % theoretical EBR of Gaussian
    theo_err_ray_prb(i)=0.5*(1-1/sqrt(1+1/SNR));     % theoretical EBR of Rayleigh

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = spectrum.welch;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%                display output 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Theoretical part

figure(1)

subplot(3,2,1);

plot(data0),title('Baseband signal');

axis([0 20000 0 0.12]);

subplot(3,2,2);

psd(h,data1,'fs',fs),title("PSD of baseband");

subplot(3,2,3);

plot(s),title('Modulated signal');

axis([0 500 -0.2 0.2]);

subplot(3,2,4);

psd(h,s,'fs',fs),title('PSD of modulated signal');

subplot(3,2,5);

plot(demodata1),title('Demodulated output');

axis([0 20000 0 1.5]);

subplot(3,2,6);

psd(h,demodata1,'fs',fs),title('PSD of demodulated output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian channel

figure(2)

subplot(2,2,1);

plot(s1),title('Modulated with Awgn');

axis([0 500 -0.5 0.5]);

subplot(2,2,2);

psd(h,s1,'fs',fs),title('PSD of modulated with Awgn)');

subplot(2,2,3);

plot(s111),title('Gaussian noise curve');

axis([0 2000 -1 1]);

subplot(2,2,4);

for i=1:nb/2

plot(idata(i),qdata(i),'ro'),title('QPSK constellation（Awgn）');hold on;

axis([-2 2 -2 2]);

plot(Awgn_ichsum(i),Awgn_qchsum(i),'b*');hold on;

legend('Theoretical（Tx）','Actual（Rx）');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gaussian followed by Rayleigh

 figure(3)

subplot(2,2,1)

plot(Ray_s),title('Modulated(Rayleigh)');

axis([0 500 -0.5 0.5]);

subplot(2,2,2);

psd(h,Ray_s,'fs',fs),title('PSD of Modulated(Ray)');

subplot(2,2,3);

for i=1:nb/2

plot(idata(i),qdata(i),'ro'),title('QPSK constellation(Rayleigh)');hold on;

axis([-2 2 -2 2]);

plot(Ray_ichsum(i),Ray_qchsum(i),'b*');hold on;

legend('Theoretical（TX）','Actual value（RX）');

end

 subplot(2,2,4)

 semilogy(SNRindB2,theo_err_awgn_prb,'r'),title('EBR');hold on;

 semilogy(SNRindB1,smld_bit_awgn_err_prb,'r*');hold on;

 semilogy(SNRindB2,theo_err_ray_prb,'b');hold on;

 semilogy(SNRindB1,smld_bit_ray_err_prb,'b*');

 xlabel('Eb/No');ylabel('BER');

 legend('theoretical AWGN','simulated AWGN','Theoretical Rayleigh','Simulated Rayleigh');

 

 



function [pb,ps]=cm_sm32(snr_in_dB)

% [pb,ps]=cm_sm32(snr_in_dB)

%                CM_SM3 finds the probability of bit error and symbol error for

%                the given value of snr_in_dB, signal to noise ratio in dB.

 

 

N=100;
E=1;                                            % energy per symbol

numofsymbolerror=0;

numofbiterror=0;

counter=0;

snr=10^(snr_in_dB/10);                          % signal to noise ratio

sgma=sqrt(E/snr)/2;                             % noise variance

s00=[1 0]; s01=[0 1]; s11=[-1 0]; s10=[0 -1];   % signal mapping

% generation of the data source

while(numofbiterror<100)

for i=1:N,

    temp=rand;                                  % a uniform random variable between 0 and 1

    if (temp<0.25),                             % with probability 1/4, source output is "00"

        dsource1(i)=0; dsource2(i)=0;

    elseif (temp<0.5),                          % with probability 1/4, source output is "01"

        dsource1(i)=0; dsource2(i)=1;

    elseif (temp<0.75),                         % with probability 1/4, source output is "10"

        dsource1(i)=1; dsource2(i)=0;

    else                                        % with probability 1/4, source output is "11"

        dsource1(i)=1; dsource2(i)=1;

    end;

end;

% detection and the probability of error calculation

 

for i=1:N,

    ray=raylrnd(0.8);

    n=sgma*randn(1,2);                          % 2 normal distributed r.v with 0, variance sgma

    if ((dsource1(i)==0) & (dsource2(i)==0)),

        r=ray*s00+n;

    elseif ((dsource1(i)==0) & (dsource2(i)==1)),

        r=ray*s01+n;

    elseif ((dsource1(i)==1) & (dsource2(i)==0)),

        r=s10*ray+n;

    else

        r=s11*ray+n;

    end;

    % The correlation metrics are computed below

    c00=dot(r,s00); c01=dot(r,s01); c10=dot(r,s10); c11=dot(r,s11);

    % The decision on the ith symbol is made next

    c_max=max([c00,c01,c10,c11]);

    if (c00==c_max), decis1=0; decis2=0;

    elseif (c01==c_max), decis1=0; decis2=1;

    elseif (c10==c_max), decis1=1; decis2=0;

    else decis1=1; decis2=1;

    end;

    % Increment the error counter, if the decision is not correct

    symbolerror=0;

    if (decis1~=dsource1(i)), numofbiterror=numofbiterror+1; symbolerror=1;

    end;

    if (decis2~=dsource2(i)), numofbiterror=numofbiterror+1; symbolerror=1;

    end;

    if (symbolerror==1), numofsymbolerror=numofsymbolerror+1;

    end;

   

end

counter=counter+1;

end

ps=numofsymbolerror/(N*counter);                          % since there are totally N symbols

pb=numofbiterror/(2*N*counter);                         % since 2N bits are transmitted

end

       

                 
function [pb1,ps1]=cm_sm33(snr_in_dB)

% [pb,ps]=cm_sm32(snr_in_dB)

%                CM_SM3 finds the probability of bit error and symbol error for

%                the given value of snr_in_dB, signal to noise ratio in dB.

N=100;
E=1;                                            % energy per symbol

snr=10^(snr_in_dB/10);                          % signal to noise ratio

sgma=sqrt(E/snr)/2;                             % noise variance

s00=[1 0]; s01=[0 1]; s11=[-1 0]; s10=[0 -1];   % signal mapping

% generation of the data source

numofsymbolerror=0;

numofbiterror=0;

counter=0;

while(numofbiterror<1000)

for i=1:N,

    temp=rand;                                  % a uniform random variable between 0 and 1

    if (temp<0.25),                             % with probability 1/4, source output is "00"

        dsource1(i)=0; dsource2(i)=0;

    elseif (temp<0.5),                          % with probability 1/4, source output is "01"

        dsource1(i)=0; dsource2(i)=1;

    elseif (temp<0.75),                         % with probability 1/4, source output is "10"

        dsource1(i)=1; dsource2(i)=0;

    else                                        % with probability 1/4, source output is "11"

        dsource1(i)=1; dsource2(i)=1;

    end;

end;

% detection and the probability of error calculation

 

for i=1:N,

    % the received signal at the detection, for the ith symbol,is:

    n=sgma*randn(1,2);                          % 2 normal distributed r.v with 0, variance sgma

    if ((dsource1(i)==0) && (dsource2(i)==0)),

        r=s00+n;

    elseif ((dsource1(i)==0) & (dsource2(i)==1)),

        r=s01+n;

    elseif ((dsource1(i)==1) & (dsource2(i)==0)),

        r=s10+n;

    else

        r=s11+n;

    end;

    % The correlation metrics are computed below

    c00=dot(r,s00); c01=dot(r,s01); c10=dot(r,s10); c11=dot(r,s11);

    % The decision on the ith symbol is made next

    c_max=max([c00,c01,c10,c11]);

    if (c00==c_max), decis1=0; decis2=0;

    elseif (c01==c_max), decis1=0; decis2=1;

    elseif (c10==c_max), decis1=1; decis2=0;

    else decis1=1; decis2=1;

    end;

    % Increment the error counter, if the decision is not correct

    symbolerror=0;

    if (decis1~=dsource1(i)), numofbiterror=numofbiterror+1; symbolerror=1;

    end;

    if (decis2~=dsource2(i)), numofbiterror=numofbiterror+1; symbolerror=1;

    end;

    if (symbolerror==1), numofsymbolerror=numofsymbolerror+1;

    end;
end
   
counter=counter+1;
end

ps1=numofsymbolerror/(N*counter);                          % since there are totally N symbols

pb1=numofbiterror/(2*N*counter);                         % since 2N bits are transmitted

end
