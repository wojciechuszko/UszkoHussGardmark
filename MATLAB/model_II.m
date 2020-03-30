%% Zooplankton community model II
% One stage-structured consumer species feeding on two resources
% For units and references, see Appendix S2, Table 2
% Created by Wojciech Uszko (2020)

function dYdt = model_II(t, Y)

%% Body masses (ng dry weight)

B_J = 100;     % juvenile
B_A = 1000;    % adult

%% Temperature- or body mass-independent parameters

deltaRS = 0.1;  % resource 1 supply rate
deltaRL = 0.1;  % resource 2 supply rate

q       = 0;    % functional response (Hill) exponent; if =0 then type II

p       = 0.85; % diet preference
pJRS    = p;
pJRL    = 1-pJRS;
pARS    = 1-pJRS;
pARL    = pJRS;

betaJ   = 0.6;  % juvenile conversion efficiency
betaA   = 0.6;  % adult conversion efficiency

HJRS    = 0.2;  % half-saturation constant
HJRL    = 0.2;  % half-saturation constant
HARS    = 0.2;  % half-saturation constant
HARL    = 0.2;  % half-saturation constant

zJA     = 0.1;   % juvenile-to-adult mass ratio

muJ     = 0.01; % juvenile background mortality rate
muA     = 0.01; % adult background mortality rate

%% Ambient temperature (Kelvin)

T = 273.15 + 20;

%% Temperature- or body mass-dependent parameters
%  Without size-temperature interaction

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = RSmax;
% 
% % Consumer maximum ingestion rate
% IJRSmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J;
% IJRLmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J;
% IARSmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;
% IARLmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;
% 
% % Consumer metabolic rate
% mJ      = (850000000 * (B_J^0.7) * exp( -0.56/(0.00008617*T) )) / B_J;
% mA      = (850000000 * (B_A^0.7) * exp( -0.56/(0.00008617*T) )) / B_A;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in temperature optimum of Imax

% Resource supply density
RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
RLmax    = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );

% Consumer maximum ingestion rate
IJRSmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_J;
IJRLmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_J;
IARSmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;
IARLmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;

% Consumer metabolic rate
mJ      = (850000000 * (B_J^0.7) * exp( -0.56/(0.00008617*T) )) / B_J;
mA      = (850000000 * (B_A^0.7) * exp( -0.56/(0.00008617*T) )) / B_A;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in metabolic rate

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );
% 
% % Consumer maximum ingestion rate
% IJRSmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J;
% IJRLmax = (19 * (B_J^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J;
% IARSmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;
% IARLmax = (19 * (B_A^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A;
% 
% % Consumer metabolic rate
% mJ      = (850000000 * (B_J^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_J;
% mA      = (850000000 * (B_A^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_A;

%% Variables

RS = Y(1);  % small resource biomass density
RL = Y(2);  % large resource biomass density
J = Y(3);   % juvenile biomass density
A = Y(4);   % adult biomass density

%% Ingestion rates

IJRS = ( pJRS * (IJRSmax/(HJRS^(1+q))) * RS^(1+q) + 0 * (IJRLmax/(HJRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJRS/(HJRS^(1+q))) * RS^(1+q) + (pJRL/(HJRL^(1+q))) * RL^(1+q) );

IJRL = ( 0 * (IJRSmax/(HJRS^(1+q))) * RS^(1+q) + pJRL * (IJRLmax/(HJRL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJRS/(HJRS^(1+q))) * RS^(1+q) + (pJRL/(HJRL^(1+q))) * RL^(1+q) );

IARS = ( pARS * (IARSmax/(HARS^(1+q))) * RS^(1+q) + 0 * (IARLmax/(HARL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pARS/(HARS^(1+q))) * RS^(1+q) + (pARL/(HARL^(1+q))) * RL^(1+q) );

IARL = ( 0 * (IARSmax/(HARS^(1+q))) * RS^(1+q) + pARL * (IARLmax/(HARL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pARS/(HARS^(1+q))) * RS^(1+q) + (pARL/(HARL^(1+q))) * RL^(1+q) );

%% Stage-specific production rates

vJ     = betaJ*(IJRS+IJRL) - mJ;                % juvenile production rate
gammaJ = (vJ - muJ) / (1 - zJA^(1-(muJ/vJ)));   % maturation rate
vA     = betaA*(IARS+IARL) - mA;                % reproduction rate

%% ODE system

dRSdt = deltaRS*(RSmax - RS) - IJRS*J - IARS*A;

dRLdt = deltaRL*(RLmax - RL) - IJRL*J - IARL*A;

dJdt = max(vA,0)*A + vJ*J - max(gammaJ,0)*J - muJ*J;

dAdt = max(gammaJ,0)*J + min(vA,0)*A - muA*A; 

dYdt = [dRSdt; dRLdt; dJdt; dAdt];

end