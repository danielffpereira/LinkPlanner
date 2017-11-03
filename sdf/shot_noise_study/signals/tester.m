clear all
close all
clc

ls=dir('*.sgn');

Power=zeros(1,length(ls));
Var=zeros(1,length(ls));

for i=1:length(ls)
    data=readSignal(ls(i).name);
    clc
    Power(i)=str2num(ls(i).name(1:length(ls(i).name)-4));
    Var(i)=var(data);
end

Power
Var

h=6.6260e-34;
c=299790000;
PdBm=-60:10;
rho=1;
g=1e6;
B=16/20e-12;
wvl=1550e-9;

VarT=1.6096261407+rho^2*g^2*B*h*c/wvl*10.^(PdBm/10)/1000;

figure
semilogy(PdBm,VarT)
hold on
semilogy(-Power,Var,'*')