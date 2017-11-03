clear all
close all
clc

color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];

ls=dir('*.sgn');

Power=zeros(1,length(ls));
Var=zeros(1,length(ls));

for i=1:length(ls)
    data=readSignal(ls(i).name);
    clc
    Power(i)=str2num(ls(i).name(1:length(ls(i).name)-4));
    Var(i)=var(data);
end

Power=-Power;


h=6.6260e-34;
c=299790000;
PdBm=-60:.001:10;
rho=1;
g=1e6;
B=16/20e-12;
wvl=1550e-9;

VarT=1.6096261407+rho^2*g^2*B*h*c/wvl*10.^(PdBm/10)/1000;

n=10.^(PdBm/10)*20e-12*wvl/(h*c)/1000;
N=10.^(Power/10)*20e-12*wvl/(h*c)/1000;

p2=1.6096261407;

%n=1:1e6;
%y=p2+p1*n;

figure
semilogy(Power,Var,'.','MarkerSize',20,'Color',color(1,:))
hold on
semilogy(PdBm,VarT,'Color',color(2,:))
semilogy([min(PdBm) max(PdBm)],[p2 p2],'--','Color',color(3,:))

xlabel('Optical Input Power [P_{dBm}]')
ylabel('Noise Variance [V^2]')
legend('Simulation Results','Theoretical Curve','Thermal Noise Level')


axis([-45 7 1 1e3])

set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')