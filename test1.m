clear;

%%%%%%%%%%%%%%%%%%%%%

%%% Initialize parameters

%%%%%%%%%%%%%%%%%%%%

T=1;                          % the width of baseband, namely the frequecy of baseband

fc=10/T;                      % carrier frequency

ml=2;                         % for separating odd and even positioning number

nb=4;                       % the number of bit to transmit

delta_T=T/100;                % sampling period

fs=1/delta_T;                 % sampling frequency

SNR=0.07;                        % signal noise ratio

t=0:delta_T:nb*T-delta_T;     % the time slot 

N=length(t);                  % the number of sampling points 
% root raised cosine parameters
rf=0.35;%roll off factor
span=5;
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

    data1((q-1)/delta_T+1:q/delta_T)=datanrz(q); % transfer polar code to waveform

end
for ii=1:N

    a(ii)=cos(2*pi*fc*t(ii));

end
%modulated signal
s=data1.*a;
figure;
plot(s);
ylim([-1.5 1.5])
title("modulated signal waveform");
%plot polar code

figure
plot(data1);
ylim([-1.5,1.5]);
title("trasmitted polar code waveform");

filtered_data=conv(data1,rcos);
figure;
plot(filtered_data);
title("filtered data waveform")

figure;
plot(abs(fft(filtered_data)));
title("filtered psd")