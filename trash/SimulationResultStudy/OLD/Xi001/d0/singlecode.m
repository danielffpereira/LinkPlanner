clear all
close all
clc

color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];

x=readSignal('MQAM4.sgn');
X=readSignal('MQAM6.sgn');

Input1=readSignal('S5.sgn');
Input2=readSignal('S6.sgn');

binary=readSignal('MQAM1.sgn');
Optical1=readSignal('MQAM8.sgn');
data1=readSignal('S9.sgn');
data2=readSignal('S12.sgn');
clc

SNU=102.5262360862843;
data1=data1/sqrt(SNU);
data2=data2/sqrt(SNU);
L=min([length(data1) length(data2)]);

state1=zeros(1,L);
state2=state1;
state3=state1;
state4=state1;

i=1;
j=1;
while i<L*2
    if binary(i)==0
        i=i+1;
        if binary(i)==0
            state1(j)=1;
        else
            state2(j)=1;
        end
    else
        i=i+1;
        if binary(i)==0
            state3(j)=1;
        else
            state4(j)=1;
        end
    end
    i=i+1;
    j=j+1;
end
state1=logical(state1);
state2=logical(state2);
state3=logical(state3);
state4=logical(state4);

power1=4*abs(Input1).^2;
power2=4*abs(Input2).^2;

load physconst
dt=20e-12/16;
wvl=1550e-9;
random1=randn(1,length(power1));
random2=randn(1,length(power2));
noise1=sqrt(h*c/(dt*wvl))*random1.*(sqrt(power1)+sqrt(h*c/(dt*wvl))*random1/4);
noise2=sqrt(h*c/(dt*wvl))*random2.*(sqrt(power2)+sqrt(h*c/(dt*wvl))*random2/4);

noisy1=power1+noise1;
noisy2=power2+noise2;
noisy=noisy1-noisy2;

clean=power1-power2;

u=1:16*15+1;


t=0:dt:(length(power1)-1)*dt;
figure
hold on

plot(t(u),clean(u),'Color',color(1,:))
plot(t(u),noisy(u),'Color',color(2,:))
xlabel('Time [s]')
ylabel('Electrical Current [A]')
legend('Without Shot Noise','With Shot Noise')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')


%%
figure
hold on
plot(data1(state1),data2(state1),'*','MarkerSize',2,'Color',color(1,:))
plot(data1(state2),data2(state2),'*','MarkerSize',2,'Color',color(2,:))
plot(data1(state3),data2(state3),'*','MarkerSize',2,'Color',color(3,:))
plot(data1(state4),data2(state4),'*','MarkerSize',2,'Color',color(4,:))
xlabel('X')
ylabel('Y')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
axis equal

figure
u=1:16*1e3;
plot(real(Optical1(u)),imag(Optical1(u)),'Color',color(1,:))

figure
u=1:16*10-1;



plot(X(u),'Color',color(1,:))
figure
plot(x(u),'Color',color(1,:))

figure
hold on
t=linspace(0,10,length(u));
plot(t,x(u),'Color',color(1,:));
plot(t,X(u),'Color',color(2,:));
legend('Before Pulse Shaper','After Pulse Shaper')
axis([0,10,-2.5,2.5])
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
xlabel('Time [T_s]')