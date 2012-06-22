function [psd] = pgen(F)
%CONSTANTS
%F=30e12;
T=5e-12;
E=10;  %The Edge parameter sets the pulse width
data=[1];
len=length(data);
%CREATE PULSE
t=-T/2:1/F:T/2;
sig=(T/8/E)^2;
if (sig==0)
    y=t;
else
    y=2*t./sig.*exp(-(t).^2/sig); 
end

%NORMALIZE PULSE AMPLITUDE TO 1
p=y./max(abs(y));

N=length(p)-1;
NEW=zeros(1,len*N+1);

ind=0;

while ind<len
    if data(ind+1)==1
        for i=(1+(ind*N)):(N+(ind*N)),
            NEW(1,i)=p(1,(i-(ind*N)));
        end
    else
        for i=(1+(ind*N)):(N+(ind*N)),
            NEW(1,i)=0;
        end
    end
    ind=ind+1;
end

%RESHAPE TIME AXIS
Tn=0:1/F:len*T;

%PLOT
plot(Tn,NEW);
grid on;

%Uncomment the following lines to obtain the PSD
%of the pulse by using Welch's method
% pwelch(NEW,128,120,length(NEW),F,'onesided')
%  psd=pwelch(NEW,128,120,length(NEW),F,'onesided');
