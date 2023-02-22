clc
clear all
close all
dp=0.0005:0.0001:0.004; %droplet diameter(m)
V=130.*(dp.^0.5); %droplet terminal velosity(m/s)
dp=dp.*100; %(cm)  
V=V.*100;   %(cm/s)
plot(dp,V,'LineWidth',3)
grid on
xlabel('droplet diameter  d_p(cm)')
ylabel('droplet terminal velocity   V_t(cm/s)')
va=0.141;   %kinematic viscisity of air at T=10 degree of celsius(cm^2/s)
Re=(V.*dp)/va;
Dg=0.141; %diffusivity of SO2 in air at T=10(cm2/s)
DL=1.38e-5; %diffusivity of SO2 in raindrop at T=10(cm2/s)
Sc=va/Dg;
if dp<0.12
    Sh=1.56+0.616.*(Re.^0.5).*(Sc^(1/3));
else
    Sh=2+0.6.*(Re.^0.5).*(Sc^(1/3));
end
sg=dp./Sh; %gas phase diffusion layer thickness(cm)
rp=dp./2; %droplet radius(cm)
sL1=0.05.*rp; %liquid phase diffusion layer thickness(cm)
sL2=0.1.*rp;
sL3=0.2.*rp;
H1=(10^(-5.6))/1000; %[H+] in normal rain(mol/cm3)
H2=(10^(-4.6))/1000; %[H+] in pre-acidified rain(mol/cm3)
kH=51.1;
k1=2.42e-5; %(mol/cm3) 
B=1/(kH*k1);
Z1=(Dg.*sL1)./(DL.*sg);
Z2=(Dg.*sL2)./(DL.*sg);
Z3=(Dg.*sL3)./(DL.*sg);
so2_g=[12 14 17 14 24 18 21 19 19 12 15 19 9 13 12 10 14 40 27 14 14 18 14 13 21 24 34 13 12 15 19 25 25 27 14 10 19 21 22 24 15 13 18 10 17 11 10 9 12 16 11 13 15 12 13 9 10 13 14 21 21 16 21 23 21 19 15 14 13 8 27 22 13 6 10 19 12 22 15 15 11 13 8]; %ppb(v)
max=max(so2_g);   mean=mean(so2_g);
so2_G=so2_g.*0.043e-12; %(mol/cm3)
SO2_g=max*0.043e-12;    %(mol/cm3) 
%for normal rain
HSO3_I11=(-1+sqrt(1+(4*B.*Z1).*(H1+Z1.*SO2_g)))./(2*B.*Z1);
HSO3_I21=(-1+sqrt(1+(4*B.*Z2).*(H1+Z2.*SO2_g)))./(2*B.*Z2);
HSO3_I31=(-1+sqrt(1+(4*B.*Z3).*(H1+Z3.*SO2_g)))./(2*B.*Z3);
FL11=(DL./sL1).*(HSO3_I11-H1); %flux in liquid layer(mole/cm2s)
FL21=(DL./sL2).*(HSO3_I21-H1);
FL31=(DL./sL3).*(HSO3_I31-H1);
figure;
subplot(1,2,1)
plot(rp,FL11,'b','LineWidth',3)
grid on
hold on 
plot(rp,FL21,'r','LineWidth',3)
plot(rp,FL31,'g','LineWidth',3)
hold off
title('normal rain:pH=5.6    SO_2=40ppb(v)')
xlabel('droplet radius:r_p(cm)');
ylabel('SO_2 flux in droplet:F_L(mole/cm^2s');
legend('s_L=5%','s_L=10%','s_L=20%')
%for pre-acidified rain
HSO3_I12=(-1+sqrt(1+(4*B.*Z1).*(H2+Z1.*SO2_g)))./(2*B.*Z1);
HSO3_I22=(-1+sqrt(1+(4*B.*Z2).*(H2+Z2.*SO2_g)))./(2*B.*Z2);
HSO3_I32=(-1+sqrt(1+(4*B.*Z3).*(H2+Z3.*SO2_g)))./(2*B.*Z3);
FL12=(DL./sL1).*(HSO3_I12-H2); %flux in liquid layer(mole/cm2s)
FL22=(DL./sL2).*(HSO3_I22-H2);
FL32=(DL./sL3).*(HSO3_I32-H2);
subplot(1,2,2)
plot(rp,FL12,'b','LineWidth',3)
grid on
hold on 
plot(rp,FL22,'r','Linewidth',3)
plot(rp,FL32,'g','LineWidth',3)
hold off
title('pre-acidified rain:pH=4.6    SO_2=40ppb(v)')
xlabel('droplet radius:r_p(cm)');
ylabel('SO_2 flux in droplet:F_L(mole/cm^2s');
legend('s_L=5%','s_L=10%','s_L=20%')