function ECON_OUT = REHEAT_ECON(ECON_IN)

power_acc = ECON_IN.power_acc; %MW, electric power provided by the accumulator alone
power_store = ECON_IN.power_store; %MW, electric power of the base case (no accumulator)
wturb = ECON_IN.wturb; %MW, turbine power output  
condenser = ECON_IN.condenser; %MW, condenser thermal power

%Price curve and amortization values
life=40; %years, amortization period
interest=0.07; %for amortization period
period=6; %hours, price period
peakAmplitude=25; %$/MWh
avgElecPrice=34; %$/MWh

%length must be connected to capital cost through pipe and insulation costs
length=ECON_IN.LTANK; %m, pipe length
Y=(1/period)*24*365; %storage cycles per year
DT = ECON_IN.discharge_time/3600; %hours, discharge time
CT = ECON_IN.charge_time/3600; %hours, charge time
n=0.8; %scaling benefits exponent, need to do linear regression

%1400 MWe Advanced Nuclear Plant Costs
ANP_turb = 809.2; %Million $
ANP_cond = 246.4; %Million $
ANP_gen = 281.4; %Million $

%Base Case (storage case)
Base_turb = ANP_turb*(power_store/1400)^n; %Million $
Base_cond = ANP_cond*(power_store/1400)^n; %Million $. 3212 MW condenser matches 1316 MW power_store
Base_gen = ANP_gen*(power_store/1400)^n; %Million $
Base_total = Base_turb+Base_cond+Base_gen; %Million $, total power dependent capital cost of base case

%Discharge case
Discharge_turb = ANP_turb*(wturb/1400)^n; %Million $
Discharge_cond = ANP_cond*(wturb/1400)^n; %Million $
Discharge_gen = ANP_gen*(wturb/1400)^n; %Million $
Discharge_total = Discharge_turb+Discharge_cond+Discharge_gen; %Million $, total power dependent capital cost of discharge case

%Charge case
Charge_turb = ANP_turb*(ECON_IN.min_load/1400)^n; %Million $
Charge_cond = ANP_cond*(ECON_IN.min_load/1400)^n; %Million $
Charge_gen = ANP_gen*(ECON_IN.min_load/1400)^n; %Million $
Charge_total = Charge_turb+Charge_cond+Charge_gen; %Million $, total power dependent capital cost of charge case

%cost of scaling up from the base case to the discharge case
Delta_cost = Discharge_total - Base_total; %Million $

%%energy dependent costs

%pipe and insulation
pipeCost=218; %$/m, cost of steel pipe
insulationCost=129.72; %$/m, cost of calcium silicate insulation. The building is insulated rather than each individual pipe.
thickness=0.203; %m, insulation thickness
pipe_and_insulation = (pipeCost*length+insulationCost*length)/10^6; %Million $

%building cost, may need to be scaled later
sqftCost=50; %$, building cost per square foot
bL=200; %m, building length
bH=10; %m, building height
bW=length/(bL*bH); %m, building width
bSA=4*(bL*bH)+2*(bW*bH); %m^2, building surface area
buildingCost=bL*bW*10.764*sqftCost/10^6; %Million $, cost to house 120 km

%storage tanks
storage_capacity = ECON_IN.storage_capacity; %m^3, required holding tank capacity
%n1 = 0.6; %use linear regression to confirm scaling factor
max_tank_capacity = 1136; %m^3
number_tanks = storage_capacity/max_tank_capacity;
storeCost = number_tanks*190.9*max_tank_capacity/10^6; %Million $

totalEnergyCost = pipe_and_insulation + buildingCost + storeCost; %Million $
totalPowerCost = Delta_cost; %Million $
totalCC=totalPowerCost+totalEnergyCost; %Million $, total overnight capital cost

c1=(3/4)*period-CT/2; %hr, charge time integral lower bound
c2=(3/4)*period+CT/2; %hr, charge time integral upper bound
d1=(period/4)-DT/2; %hr, discharge time integral lower bound
d2=(period/4)+DT/2; %hr, discharge time integral upper bound
y=@(t)peakAmplitude*sin((2*pi()*t)/period)+avgElecPrice; 
intC=integral(y,c1,c2); %$/MW
intD=integral(y,d1,d2); %$/MW
ADP= intD/(d2-d1); %$/MWh, Average discharge price
ACP= intC/(c2-c1); %$/MWh, Average charge price
DP=ADP-ACP ; %$/MWh, delta price

RC=ACP*CT*Y*(power_store-ECON_IN.min_load)/10^6; %M$/year, forgone revenue from charging
RD=ADP*DT*Y*power_acc/10^6; %M$/year, revenue from discharging
CC=totalCC*(interest+(interest/((1+interest)^life-1))); %M$/year, amortized capital cost
totalOM = 0.05*CC; %M$/year, O&M estimated as 5% yearly capital cost
netRevenue=RD-RC-CC-totalOM; %M$/year, revenue provided by the addition of the accumulator

ECON_OUT.netRevenue = netRevenue; %M$/year
ECON_OUT.CC = CC; %M$/year
ECON_OUT.RC = RC; %M$/year
ECON_OUT.RD = RD; %M$/year
ECON_OUT.OM = totalOM; %M$/year
ECON_OUT.totalCC = totalCC; %M$, total overnight capital cost
%ECON_OUT.revenueEff = revenueEff;

end



