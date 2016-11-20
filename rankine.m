function outdata = rankine(indata)
%RANKINE A simulation of the ideal Rankine Cycle
%   RANKINE(IN) generates the thermodynamic properties of the ideal
%   rankine cycle.  
%
%   Input parameters:
%        in.p1    - High side pressure on the turbine (bar)
%        in.p2    - Low  side pressure on the turbine (bar)
%        in.Wd    - Turbine output power in W
%        in.t2    - Coolant  inlet temperature
%        in.t3    - Coolant outlet temperature
%
%   Output parameters:
%
%       out.md    - mass flow rate (kg/s)
%       out.mu    - thermodynamic efficiency
%       out.bwr   - back work ratio
%       out.Qdin  - Rate of energy in (W)
%       out.Qdout - Rate of energy out (W)
%       out.mdcw  - Condenser mass flow rate
%
% This code requires XSteam.m, which available from MATLAB Central.

% Michael Agostini
% Copyright 2006-2009 The MathWorks, Inc.
% This is provided as an example only.


% Load the data from a structure into individual variable names

   p1  = indata.p1;
   p2  = indata.p2;
   Wd  = indata.Wd;
   t2  = indata.t2;
   t3  = indata.t3;

%% Look up thermodynamic properties at State 1

   h1  = XSteam('hV_p', p1 )  ;
   s1  = XSteam('sV_p', p1 )  ;

%% Look up and compute properties for State 2
   
   s2  =  s1                  ;

   sg  = XSteam('sV_p', p2 )  ;
   sf  = XSteam('sL_p', p2 )  ;

   hg  = XSteam('hV_p', p2 )  ;
   hf  = XSteam('hL_p', p2 )  ;
   hfg = hg - hf              ;

   x2  = (s2-sf) / (sg-sf)    ;

   h2  = hf + x2*hfg          ;

%% State 3 properties

   h3  = hf                   ;

%% State 4 properties

   v3  = XSteam('vL_p', p2 )  ;
   h4  = h3 + v3*(p1 - p2)*1e3;

%% Compute the Thermodynamic Efficiency 

   mu  = ((h1-h2) - (h4-h3)) / (h1-h4);

%% Compute the Backwork Ratio

   bwr = (h4-h3) / (h1-h2) ;

%% Compute the Mass flow rate at the condenser

   md  =  Wd/( ((h1-h2)-(h4-h3))*1000 ) ;

%% Qd (Energy flow)

   Qdin  = md * (h1-h4)*1000 ; 
   Qdout = md * (h2-h3)*1000 ;

%% Steady state energy

   hC    = XSteam('hL_T', t2      )  ;
   hH    = XSteam('hL_T', t3      )  ;

   mdcw  = md * (h2-h3)/(hH-hC)      ;

%% Pack the Output data into a structure

   outdata.md    = md           ;
   outdata.mu    = mu           ;
   outdata.bwr   = bwr          ;
   outdata.Qdin  = Qdin         ;
   outdata.Qdout = Qdout        ;
   outdata.mdcw  = mdcw         ;
   
 
