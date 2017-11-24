clear all
close all
clc

ls=dir;
ls(1)=[];
ls(1)=[];


binary=readSignal('MQAM1.sgn');
data1=readSignal('S9.sgn');
data2=readSignal('S12.sgn');
clc

L=min([length(data1) length(data2)]);
data1=data1(1:L);
data2=data2(1:L);

SNU=102.5262360862843;
n=.5*2;
VA=.5*2;

y1=data1/(sqrt(SNU));
y2=data2/(sqrt(SNU));

x1=zeros(1,L);
x2=x1;
state1=x1;
state2=x1;
state3=x1;
state4=x1;

i=1;
j=1;

D=1;
while i<L*2
    if binary(i)==0
        i=i+1;
        if binary(i)==0
            x1(j)= sqrt(n);
            x2(j)=-sqrt(n);
            state1(j)=1;
        else
            x1(j)=-sqrt(n);
            x2(j)= sqrt(n);
            state2(j)=1;
        end
    else
        i=i+1;
        if binary(i)==0
            x1(j)=-sqrt(n);
            x2(j)=-sqrt(n);
            state3(j)=1;
        else
            x1(j)= sqrt(n);
            x2(j)= sqrt(n);
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

t1=sum(y1.*x1)/sum(x1.^2);
t2=sum(y2.*x2)/sum(x2.^2);
t=mean([t1 t2]);
T=t1^2+t2^2;

sigma1=1/L*sum((y1-t1*x1).^2)*4;
sigma2=1/L*sum((y2-t2*x2).^2)*4;
sigma=mean([sigma1 sigma2]);


epe=1e-10;
zepe=sqrt(2)*erfinv(1-epe);
tmin=t-zepe*sqrt(sigma/(L*VA));
sigmamax=sigma*(1+zepe*sqrt(2)/sqrt(L));

lambda=[1/2*exp(-VA/2)*(cosh(VA/2)+cos(VA/2)) ...
    1/2*exp(-VA/2)*(sinh(VA/2)+sin(VA/2)) ...
    1/2*exp(-VA/2)*(cosh(VA/2)-cos(VA/2)) ...
    1/2*exp(-VA/2)*(sinh(VA/2)-sin(VA/2))];

Z=VA*(lambda(1)^(3/2)/lambda(2)^(1/2)+...
    lambda(2)^(3/2)/lambda(3)^(1/2)+...
    lambda(3)^(3/2)/lambda(4)^(1/2)+...
    lambda(4)^(3/2)/lambda(1)^(1/2));

Id=eye(2);
g=[1 0;0 -1];

SNR=t^2*VA/(sigma);
I=log2(1+SNR);

Gamma1=(VA+1)*Id;
Gamma2=(t^2*VA+sigma)*Id;
C12=t*Z*g;
Gamma=[Gamma1 C12; C12 Gamma2 ];
GammaAB=(VA+1-t^2*Z^2/(t^2*VA+1+sigma))*Id;

Delta=det(Gamma1)+det(Gamma2)+2*det(C12);
k1=eig(GammaAB);
mu1=sqrt(.5*(Delta+sqrt(Delta^2-4*det(Gamma))));
mu2=sqrt(.5*(Delta-sqrt(Delta^2-4*det(Gamma))));
mu3=sqrt(k1(1)*k1(2));
g= @(x) ((x+1)/2).*log2((x+1)/2)-((x-1)/2).*log2((x-1)/2);
S=g(mu1)+g(mu2)-g(mu3);

beta=.8;
K=(beta*I-S)
XXX=(sigma-1)*2/T


%%
clear all
close all
clc

color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];

data1=readSignal('S9.sgn');
data2=readSignal('S12.sgn');
binary=readSignal('MQAM1.sgn');
clc
L=min([length(data1) length(data2)]);

SNU=102.5262360862843;
VA=.5*2;

data1=data1/sqrt(SNU);
data2=data2/sqrt(SNU);

x1=zeros(1,L);
x2=x1;
state1=x1;
state2=x1;
state3=x1;
state4=x1;

i=1;
j=1;

D=1;
while i<L*2
    if binary(i)==0
        i=i+1;
        if binary(i)==0
            x1(j)= sqrt(VA/D);
            x2(j)=-sqrt(VA/D);
            state1(j)=1;
        else
            x1(j)=-sqrt(VA/D);
            x2(j)= sqrt(VA/D);
            state2(j)=1;
        end
    else
        i=i+1;
        if binary(i)==0
            x1(j)=-sqrt(VA/D);
            x2(j)=-sqrt(VA/D);
            state3(j)=1;
        else
            x1(j)= sqrt(VA/D);
            x2(j)= sqrt(VA/D);
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

figure
hold on
plot(data1(state1),data2(state1),'*','MarkerSize',2,'Color',color(1,:))
plot(data1(state2),data2(state2),'*','MarkerSize',2,'Color',color(2,:))
plot(data1(state3),data2(state3),'*','MarkerSize',2,'Color',color(3,:))
plot(data1(state4),data2(state4),'*','MarkerSize',2,'Color',color(4,:))
xlabel('I [SNU]')
ylabel('Q [SNU]')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
axis equal
axis([-1 1 -1 1]*3)