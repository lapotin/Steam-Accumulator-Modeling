function outdata = RANKINE_REHEAT(in)

ELECTRICPOWERACC = in.Power_Electric; %MW, desired electric power of the accumulator
ELECTRICPOWERSTORECASE = 1316; %MW, electrical power in storage case (rankine reheat without accumulator contribution)
MDOTSTORE = 1333; %kg/s, mass flow rate w/ out accumulator
ELECTRICPOWER = ELECTRICPOWERACC + ELECTRICPOWERSTORECASE; %MW, rankine cycle total electric power

%outlet of steam generator (inlet to HPT)
p2 = 72; %bar
t2 = 288; %C
h2 = 2770; %kJ/kg
x2 = 1;

THERMALPOWER = 3500; %MW, constant to maintain reactor
APPARENTEFF = ELECTRICPOWER/THERMALPOWER;

MDOT = (ELECTRICPOWER/ELECTRICPOWERSTORECASE)*MDOTSTORE; %kg/s, mass flow rate of feedwater heater

%inlet to steam generator
p1 = 72; %bar
h1 = h2-(THERMALPOWER*1000)/MDOT; %kJ/kg
t1 = XSteam('T_ph',p1,h1); %C

%outlet of pump (inlet to feedwater heater)
p6=72; %bar
t6=39; %C
h6=171; %kJ/kg

QFEEDWATERHEATER = MDOT*(h1-h6)/1000; %MW
%outlet of steam accumulator discharge after feedwater heater
p7 = 1; %bar
t7 = 49; %C
h7 = 206; %kJ/kg
rho7 = XSteam('rho_pT',p7,t7); %kg/m^3
hACC = in.hACC; %kJ/kg
MDOTACC = (QFEEDWATERHEATER*1000)/(hACC-h7); %kg/s, discharge rate of accumulator

%outlet of condenser (inlet to pump)
x5 = 0; %inlet to pump is a saturated liquid
t5 = 39; %C
p5 = XSteam('psat_T',t5); %bar
h5 = XSteam('hL_p',p5); %kJ/kg
PUMP=MDOT*(h6-h5)/1000; %MW

WTURB = ELECTRICPOWER + PUMP; %MW, total turbine work

%outlet of LPT
h4 = h2-(WTURB*1000)/MDOT; %kJ/kg
CONDENSER = MDOT*(h4-h5)/1000; %MW

ACCEFF = ELECTRICPOWERACC/QFEEDWATERHEATER;
PERCENTLOAD = (WTURB/(ELECTRICPOWERSTORECASE+PUMP))*100; %percent of base electric output by turbine
EFF1 = ELECTRICPOWER/(THERMALPOWER + QFEEDWATERHEATER);

outdata.MDOTACC = MDOTACC;
outdata.MDOT = MDOT;
outdata.QFEEDWATERHEATER = QFEEDWATERHEATER;
outdata.PUMP = PUMP;
outdata.ELECTRICPOWERSTORECASE = ELECTRICPOWERSTORECASE;
outdata.MDOTSTORE = MDOTSTORE;
outdata.THERMALPOWER = THERMALPOWER;
outdata.ELECTRICPOWER = ELECTRICPOWER;
outdata.APPARENTEFF = APPARENTEFF;
outdata.ACCEFF = ACCEFF;
outdata.WTURB = WTURB;
outdata.PERCENTLOAD = PERCENTLOAD;
outdata.CONDENSER = CONDENSER;
outdata.h4 = h4; 
outdata.h6 = h6;
outdata.rho7 = rho7;
outdata.EFF1= EFF1;

end
