clear all;
close all;

PASCALS_PER_BAR = 1.e5;
EPSILON = .01;

DT = 10;  % s
P0 = 70; % bar
Power = 500; %MWe, turbine power output
Energy = 1500; %MWh, turbine energy
RTANK = 0.4064; % m (16 inches)
LTANK = 338000.; % m
VTANK = LTANK.*pi.*RTANK^2.;  % m3
T_END = (Energy/Power)*3600; % run time (s)
X0 = 0.06; % vapor quality (mass fraction) 
TAU = 85; % s, relaxation time
deltaH=-.01; %kJ/kcg, small change in enthalpy used to evaluate numerical derivative
deltaP=-.01; %bar, small change in pressure used to evaluate numerical derivative

%QLOSS=.286*LTANK; %kW, heat loss calculated in Qloss to be 286 W/m
QLOSS=0;

%mass flow rate of exiting vapor is calculated in the Rankine calculator further down in
%this model
VAP_IN = 0; 
%VAP_OUT = ;
LIQ_IN = 0;
LIQ_OUT = 0;

PVAP_IN = 70.; %bar
PLIQ_IN = 70.;
HVAP_IN = XSteam('hV_p',PVAP_IN); %kj/kg
HLIQ_IN = XSteam('hL_p',PLIQ_IN); %kj/kg

N = round(T_END/DT)+1;
time=zeros(N,1);
time(1)=0.;
for(i=1:N-1)
    time(i+1)=time(i)+DT;
end

p=zeros(N,1);  % bar
m1=zeros(N,1);  % kg
m2=zeros(N,1);
Vol1=zeros(N,1);  % m3
Vol2=zeros(N,1); %m3
v1=zeros(N,1); %m3/kg
v2=zeros(N,1); %m3/kg
rho1=zeros(N,1);   % kg/m3
rho2=zeros(N,1);
h1=zeros(N,1);  % kJ/kg 
h2=zeros(N,1);
h2in=zeros(N,1);
h2out=zeros(N,1);
x=zeros(N,1);
x_1=zeros(N,1); %x1 test
x_2=zeros(N,1); %x2 test
t1=zeros(N,1);   % C
t2=zeros(N,1);
m1b=zeros(N,1);   % kg/s
m2b=zeros(N,1);
mh1b=zeros(N,1);  % kJ/s
mh2b=zeros(N,1);
mpt1=zeros(N,1);   % kg/s
mpt2=zeros(N,1);   % kg/s; mpt2=-mpt1

total_m2_discharged=0;
average_efficiency=0;
% initial setup. 
i=1;
x(i)=X0; %initial quality in tank
x_1(i)=0;
x_2(i)=1;
p(i)=P0; %initial pressure in tank (bar)
rho1(i)=XSteam('rhoL_p',p(i)); %kg/m3
rho2(i)=XSteam('rhoV_p',p(i)); %kg/m3
t1(i)=XSteam('Tsat_p',p(i));
t2(i)=t1(i);
mix_rho=1./(x(i)/rho2(i)+(1.-x(i))/rho1(i));
m_total=mix_rho*VTANK; %total mass in tank
m1(i)=m_total*(1-x(i)); %total liquid mass in tank (kg)
m2(i)=m_total*x(i); %total vapor mass in tank (kg)
Vol1(i)=m1(i)/rho1(i); %liquid volume in tank (m3)
Vol2(i)=m2(i)/rho2(i); %vapor volume in tank (m3)
Vol_total(i)=Vol1(i)+Vol2(i); %volume of liquid and vapor in tank should be conserved(m3)
v1(i)=1/rho1(i); %liquid specific volume in tank (m3/kg)
v2(i)=1/rho2(i); %vapor specific volume in tank (m3/kg)
h1(i) = XSteam('h_pT',p(i),t1(i)-EPSILON); %(kj/kg)
h2(i) = XSteam('h_pT',p(i),t2(i)+EPSILON);

Dh=(XSteam('v_ph',p(i),h1(i)+deltaH)-XSteam('v_ph',p(i),h1(i)))/abs(deltaH); %numerical derivative
Dp=(XSteam('v_ph',p(i)+deltaP,h1(i))-XSteam('v_ph',p(i),h1(i)))/abs(deltaP); %numerical derivative

Qloss1(i)=(Vol1(i)/(Vol1(i)+Vol2(i)))*QLOSS; %volume fraction of heat lost by liquid (kW)
Qloss2(i)=(Vol2(i)/(Vol1(i)+Vol2(i)))*QLOSS; %volume fraction of heat lost by vapor (kW)

%Rankine properties

in.p1=P0; %high turbine pressure (bar)
in.p2=0.9; %low turbine pressure (bar)
in.Wd=Power*10^6; %turbine output power (W)
in.t2=15;
in.t3=40;
outdata(i)=rankine(in);
eff(i)=outdata(i).mu;
VAP_OUT(i)=outdata(i).md; %mass flow rate of exiting vapor



for i=1:(N-1)
    % calculate initial state

    m1b(i)=LIQ_IN-LIQ_OUT; %kg/s
    m2b(i)=VAP_IN-VAP_OUT(i); %kg/s
    total_m2_discharged=total_m2_discharged+VAP_OUT(i)*DT;
    
    mh1b(i)=LIQ_IN*HLIQ_IN-LIQ_OUT*h1(i); %kW
    mh2b(i)=VAP_IN*HVAP_IN-VAP_OUT(i)*h2(i); %kW

    hL_sat = XSteam('hL_p',p(i)); %saturated liquid enthalpy (kJ/kg)
    hV_sat = XSteam('hV_p',p(i)); %saturated vapor enthalpy (kJ/kg)

    r = hV_sat - hL_sat;     % latent heat of vaporization
    if(h1(i) > hL_sat)
        mc = 0.;
        me(i) = (rho1(i)*Vol1(i)*(h1(i)-hL_sat))/(TAU*r);
    else
        me(i) = 0.;
        mc = rho1(i)*Vol1(i)*(hL_sat-h1(i))/(TAU*r);
    end
    
    mpt1(i)=mc-me(i); %kg/s, liquid mass change due to evaporation and condensation
    mpt2(i)=me(i)-mc; %kg/s, vapor mass change due to evaporation and condensation
    m1(i+1)=m1(i)+(m1b(i)+mpt1(i))*DT; %liquid mass balance (kg)
    m2(i+1)=m2(i)+(m2b(i)+mpt2(i))*DT; %vapor mass balance (kg)
    
    dv1dh=(XSteam('v_ph',p(i),h1(i)+Dh)-XSteam('v_ph',p(i),h1(i)))/Dh; % (m3/kJ), change in v1 per change in h at constant p
    dv2dh=(XSteam('v_ph',p(i),h2(i)+Dh)-XSteam('v_ph',p(i),h2(i)))/Dh; % (m3/kJ)
    dv1dp=(XSteam('v_ph',p(i)+Dp,h1(i))-XSteam('v_ph',p(i),h1(i)))/(Dp*PASCALS_PER_BAR); %(m5/(N*kg)) change in v1 per change in p at constant h
    dv2dp=(XSteam('v_ph',p(i)+Dp,h2(i))-XSteam('v_ph',p(i),h2(i)))/(Dp*PASCALS_PER_BAR);
    
    term1 = (h1(i)*dv1dh-v1(i))*(m1(i+1)-m1(i))/DT; %m3/s
    term2 = (h2(i)*dv2dh-v2(i))*(m2(i+1)-m2(i))/DT; %m3/s
    term3 = dv1dh*(mh1b(i)+mpt1(i)*hV_sat-Qloss1(i)); %m3/s
    term4 = dv2dh*(mh2b(i)+mpt2(i)*hV_sat-Qloss2(i)); %m3/s
    term5 = m1(i)*(dv1dp+v1(i)*dv1dh*.001); %m5/N
    term6 = m2(i)*(dv2dp+v2(i)*dv2dh*.001); %m5/N
    
    p(i+1) = ((DT*((term1+term2-term3-term4)/(term5+term6)))/PASCALS_PER_BAR)+p(i); %bar
    h1(i+1) = (mh1b(i)+mpt1(i)*hV_sat+(Vol1(i)/DT)*PASCALS_PER_BAR*(p(i+1)-p(i))*.001-(h1(i)/DT)*(m1(i+1)-m1(i)))*(DT/m1(i))+h1(i); %kJ/kg
    h2(i+1) = (mh2b(i)+mpt2(i)*hV_sat+(Vol2(i)/DT)*PASCALS_PER_BAR*(p(i+1)-p(i))*.001-(h2(i)/DT)*(m2(i+1)-m2(i)))*(DT/m2(i))+h2(i); %kJ/kg
    t1(i+1) = XSteam('T_ph',p(i+1),h1(i+1));
    t2(i+1) = XSteam('T_ph',p(i+1),h2(i+1));
    rho1(i+1) = XSteam('rho_ph',p(i+1),h1(i+1));
    rho2(i+1) = XSteam('rho_ph',p(i+1),h2(i+1));
    Vol1(i+1) = m1(i+1)/rho1(i+1);
    Vol2(i+1) = m2(i+1)/rho2(i+1); 
    Vol_total(i+1)=Vol1(i+1)+Vol2(i+1); 
    v1(i+1) = 1/rho1(i+1);
    v2(i+1) = 1/rho2(i+1);
    x(i+1) = m2(i+1)/(m1(i+1)+m2(i+1));
    x_1(i+1) = XSteam('vx_ph',p(i+1),h1(i+1));
    x_2(i+1) = XSteam('vx_ph',p(i+1),h2(i+1));
    Qloss1(i+1)=(Vol1(i+1)/(Vol1(i+1)+Vol2(i+1)))*QLOSS; %(kW)
    Qloss2(i+1)=(Vol2(i+1)/(Vol1(i+1)+Vol2(i+1)))*QLOSS; %(kW)
    
    
    %New Rankine Properties
    in.p1=p(i+1);
    outdata(i+1)=rankine(in);
    eff(i+1)=outdata(i+1).mu;
    VAP_OUT(i+1)=outdata(i+1).md;
    average_efficiency=average_efficiency+DT*eff(i+1);
end

percentVolError=((Vol_total(i+1)-Vol_total(1))/Vol_total(i+1))*100;
%average_efficiency=average_efficiency*(1/T_END);

% figure(1)
% plot(time,p);
% hold on;
% xlabel('Time [s]');
% ylabel('Pressure [bar]');
% 
% figure(2)
% plot(time,m1,'r');
% hold on;
% plot(time,m2);
% xlabel('Time [s]');
% ylabel('Mass in accumulator [kg]');
% legend('Liquid','Vapor');
% 
% figure(3)
% plot(time,t1,'r');
% hold on;
% plot(time,t2);
% xlabel('Time [s]');
% ylabel('Temperatue [C]');
% legend('Liquid','Vapor');
% 
% figure(4)
% plot(time,Vol1,'r');
% hold on;
% plot(time,Vol2);
% plot(time,Vol1+Vol2,'g');
% xlabel('Time [s]');
% ylabel('Phase Volume [m3]');
% legend('Liquid','Vapor','Total');
% 
% figure(5)
% plot(time,eff*100);
% hold on;
% xlabel('Time [s]');
% ylabel('Efficiency [%]');
% 
% figure(6)
% plot(time,VAP_OUT);
% hold on;
% xlabel('Time [s]');
% ylabel('Mass Flow Rate [kg/s]');

figure(7)
plot(time, Vol_total,'g');
hold on;
plot(time,VTANK,'r');
xlabel('Time [s]');
ylabel('Volume [m3]');
legend('phase volumes','volume of tank');

figure(8)
plot(time,x*100,'g');
hold on;
plot(time,x_1*100,'r');
hold on;
plot(time,x_2*100,'b');
xlabel('Time [s]');
ylabel('Quality [%]');
legend('total quality','quality of liquid','quality of vapor');

