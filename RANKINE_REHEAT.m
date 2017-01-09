function outdata = RANKINE_REHEAT(in)

ELECTRICPOWERACC = in.Power_Electric; %MW, desired electric power of the accumulator
MDOT = 1744; %kg/s, mass flow rate of feedwater heater
ELECTRICPOWERSTORECASE = 1316; %MW, electrical power in storage case (rankine reheat without accumulator contribution)
MDOTSTORE = 1333; %kg/s, mass flow rate w/ out accumulator
ELECTRICPOWER = ELECTRICPOWERACC + ELECTRICPOWERSTORECASE; %MW, rankine cycle total electric power

%inlet to steam generator
p1 = 72; %bar
t1 = 179; %C
h1 = 763; %kJ/kg

%outlet of steam generator (inlet to HPT)
p2 = 72; %bar
t2 = 288; %C
h2 = 2770; %kJ/kg
x2 = 1;

THERMALPOWER = MDOT*(h2-h1)/1000; %MW
APPARENTEFF = ELECTRICPOWER/THERMALPOWER;

%assume ELECTRICPOWER ~ WTURB
WTURB = ELECTRICPOWER; %MW, total turbine work

%outlet of LPT
h4 = h2-(WTURB*1000)/MDOT; %kJ/kg

%outlet of condenser (inlet to pump)
x5 = 0; %inlet to pump is a saturated liquid
t5 = 39; %C
p5 = XSteam('psat_T',t5); %bar
h5 = XSteam('hL_p',p5); %kJ/kg
CONDENSER = MDOT*(h4-h5)/1000; %MW
QFEEDWATERHEATER = ELECTRICPOWER-THERMALPOWER+CONDENSER; %MW

%outlet of pump (inlet to feedwater heater)
h6 = h1-(QFEEDWATERHEATER*1000)/MDOT; %kJ/kg

%outlet of steam accumulator discharge after feedwater heater
p7 = 1; %bar
t7 = 49; %C
h7 = 206; %kJ/kg
rho7 = XSteam('rho_pT',p7,t7); %kg/m^3

hACC = in.hACC; %kJ/kg
MDOTACC = (QFEEDWATERHEATER*1000)/(hACC-h7); %kg/s, discharge rate of accumulator

ACCEFF = ELECTRICPOWERACC/QFEEDWATERHEATER;
PERCENTLOAD = ELECTRICPOWERSTORECASE/ELECTRICPOWER; %percent of based electric output
EFF1 = ELECTRICPOWER/(THERMALPOWER + QFEEDWATERHEATER);

outdata.MDOTACC = MDOTACC; 
outdata.QFEEDWATERHEATER = QFEEDWATERHEATER;
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
