clear;
 close all;
 clc

load('figure5a.mat')

plotimage(seis_record)
seis_recorda=seis_record;


load('figure5b.mat')
plotimage(seis_record)
seis_recordb=seis_record;



load('figure5c.mat')
plotimage(    seis_record)
seis_recordc=seis_record;


load('figure5d.mat')
plotimage(    seis_record)
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,50);


figure;plot(seis_recorda(:,50),'b','linewidth',2);

hold on;plot(seis_recordb(:,50),'r','linewidth',2)
hold on;plot(seis_recordc(:,50),'k','linewidth',2)

hold on;plot(temp,'m','linewidth',2)

% hold on;plot(temp-1*10^-5,'b','linewidth',1);

xlabel('time(ms)')
ylabel('Amp')
legend('FD scheme without L-W','FD scheme with L-W','proposed FD scheme with L-W','Time dispersion elimination')
grid on
% axis([ 0 1500 -7*10^-5 10*10^-5])
