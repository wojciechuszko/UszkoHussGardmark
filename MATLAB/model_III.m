%% Zooplankton community model III
% Two stage-structured consumer species feeding on two resources
% For units and references, see Appendix S2, Table 2
% Created by Wojciech Uszko (2020)

function dYdt = model_III(t, Y)

%% Body masses (ng dry weight)

B_J1 = 100;     % juvenile of consumer 1
B_A1 = 1000;    % adult of consumer 1
B_J2 = 1000;    % juvenile of consumer 2
B_A2 = 10000;   % adult of consumer 1

%% Temperature- or body mass-independent parameters

deltaRS = 0.1;  % resource 1 supply rate
deltaRL = 0.1;  % resource 2 supply rate

q       = 0;    % functional response (Hill) exponent; if =0 then type II

p       = 0.75; % diet preference
pA1RS   = p;
pA1RL   = 1-pA1RS;
pJ2RL   = pA1RS;
pJ2RS   = 1-pA1RS;

betaJ1   = 0.6;  % juvenile 1 conversion efficiency
betaA1   = 0.6;  % adult 1 conversion efficiency
betaJ2   = 0.6;  % juvenile 2 conversion efficiency
betaA2   = 0.6;  % adult 2 conversion efficiency

HJ1RS   = 0.2;  % half-saturation constant
HA1RS   = 0.2;  % half-saturation constant
HA1RL   = 0.2;  % half-saturation constant
HJ2RS   = 0.2;  % half-saturation constant
HJ2RL   = 0.2;  % half-saturation constant
HA2RL   = 0.2;  % half-saturation constant

zJ1A1   = 0.1;  % juvenile-to-adult mass ratio in consumer 1
zJ2A2   = 0.1;  % juvenile-to-adult mass ratio in consumer 2

muJ1    = 0.01; % juvenile 1 background mortality rate
muA1    = 0.01; % adult 1 background mortality rate
muJ2    = 0.01; % juvenile 2 background mortality rate
muA2    = 0.01; % adult 2 background mortality rate

%% Ambient temperature (Kelvin)

T = 273.15 + 20;

%% Temperature- or body mass-dependent parameters
%  Without size-temperature interaction

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = RSmax;
% 
% % Consumer maximum ingestion rate
% IJ1RSmax = (19 * (B_J1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J1;
% IA1RSmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
% IA1RLmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
% IJ2RSmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
% IJ2RLmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
% IA2RLmax = (19 * (B_A2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A2;
% 
% % Consumer metabolic rate
% mJ1      = (850000000 * (B_J1^0.7) * exp( -0.56/(0.00008617*T) )) / B_J1;
% mA1      = (850000000 * (B_A1^0.7) * exp( -0.56/(0.00008617*T) )) / B_A1;
% mJ2      = (850000000 * (B_J2^0.7) * exp( -0.56/(0.00008617*T) )) / B_J2;
% mA2      = (850000000 * (B_A2^0.7) * exp( -0.56/(0.00008617*T) )) / B_A2;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in temperature optimum of Imax

% Resource supply density
RSmax    = 0.0042 * exp( 0.151/(0.00008617*T) );
RLmax    = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );

% Consumer maximum ingestion rate
IJ1RSmax = (19 * (B_J1^(0.7)) * exp(-((T-(273.15+24))^2)/(2*(8^2)))) / B_J1;
IA1RSmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
IA1RLmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
IJ2RSmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
IJ2RLmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
IA2RLmax = (19 * (B_A2^(0.7)) * exp(-((T-(273.15+16))^2)/(2*(8^2)))) / B_A2;

% Consumer metabolic rate
mJ1      = (850000000 * (B_J1^0.7) * exp( -0.56/(0.00008617*T) )) / B_J1;
mA1      = (850000000 * (B_A1^0.7) * exp( -0.56/(0.00008617*T) )) / B_A1;
mJ2      = (850000000 * (B_J2^0.7) * exp( -0.56/(0.00008617*T) )) / B_J2;
mA2      = (850000000 * (B_A2^0.7) * exp( -0.56/(0.00008617*T) )) / B_A2;

%% Temperature- or body mass-dependent parameters
%  With size-temperature interaction in Rmax and in metabolic rate

% % Resource supply density
% RSmax   = 0.0042 * exp( 0.151/(0.00008617*T) );
% RLmax   = (5.88* 10^(-7)) * exp( 0.37564/(0.00008617*T) );
% 
% % Consumer maximum ingestion rate
% IJ1RSmax = (19 * (B_J1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J1;
% IA1RSmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
% IA1RLmax = (19 * (B_A1^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A1;
% IJ2RSmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
% IJ2RLmax = (19 * (B_J2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_J2;
% IA2RLmax = (19 * (B_A2^(0.7)) * exp(-((T-(273.15+20))^2)/(2*(8^2)))) / B_A2;
% 
% % Consumer metabolic rate
% mJ1      = (850000000 * (B_J1^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_J1;
% mA1      = (850000000 * (B_A1^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_A1;
% mJ2      = (850000000 * (B_J2^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_J2;
% mA2      = (850000000 * (B_A2^(0.7 + 0.0005*T)) * exp( -0.56/(0.00008617*T) )) / B_A2;

%% Variables

RS = Y(1);  % small resource biomass density
RL = Y(2);  % large resource biomass density
J1 = Y(3);  % juvenile 1 biomass density
A1 = Y(4);  % adult 1 biomass density
J2 = Y(5);  % juvenile 2 biomass density 
A2 = Y(6);  % adult 1 biomass density

%% Ingestion rates

IJ1RS = ( 1 * (IJ1RSmax/(HJ1RS^(1+q))) * RS^(1+q) ) / ... 
    ( 1 + (1/(HJ1RS^(1+q))) * RS^(1+q) );

IA1RS = ( pA1RS * (IA1RSmax/(HA1RS^(1+q))) * RS^(1+q) + 0 * (IA1RLmax/(HA1RL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pA1RS/(HA1RS^(1+q))) * RS^(1+q) + (pA1RL/(HA1RL^(1+q))) * RL^(1+q) );

IA1RL = ( 0 * (IA1RSmax/(HA1RS^(1+q))) * RS^(1+q) + pA1RL * (IA1RLmax/(HA1RL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pA1RS/(HA1RS^(1+q))) * RS^(1+q) + (pA1RL/(HA1RL^(1+q))) * RL^(1+q) );

IJ2RS = ( pJ2RS * (IJ2RSmax/(HJ2RS^(1+q))) * RS^(1+q) + 0 * (IJ2RLmax/(HJ2RL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJ2RS/(HJ2RS^(1+q))) * RS^(1+q) + (pJ2RL/(HJ2RL^(1+q))) * RL^(1+q) );

IJ2RL = ( 0 * (IJ2RSmax/(HJ2RS^(1+q))) * RS^(1+q) + pJ2RL * (IJ2RLmax/(HJ2RL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (pJ2RS/(HJ2RS^(1+q))) * RS^(1+q) + (pJ2RL/(HJ2RL^(1+q))) * RL^(1+q) );

IA2RL = ( 1 * (IA2RLmax/(HA2RL^(1+q))) * RL^(1+q) ) / ... 
    ( 1 + (1/(HA2RL^(1+q))) * RL^(1+q) );

%% Stage-specific production rates

vJ1     = betaJ1*IJ1RS - mJ1;                           % juvenile 1 production rate
gammaJ1 = (vJ1 - muJ1) / (1 - zJ1A1^(1-(muJ1/vJ1)));    % maturation rate
vA1     = betaA1*(IA1RS+IA1RL) - mA1;                   % reproduction rate

vJ2     = betaJ2*(IJ2RS+IJ2RL) - mJ2;                   % juvenile 2 production rate
gammaJ2 = (vJ2 - muJ2) / (1 - zJ2A2^(1-(muJ2/vJ2)));    % maturation rate
vA2     = betaA2*IA2RL - mA2;                           % reproduction rate

%% ODE system

dRSdt = deltaRS*(RSmax - RS) - IJ1RS*J1 - IA1RS*A1 - IJ2RS*J2;

dRLdt = deltaRL*(RLmax - RL) - IA1RL*A1 - IJ2RL*J2 - IA2RL*A2;

dJ1dt = max(vA1,0)*A1 + vJ1*J1 - max(gammaJ1,0)*J1 - muJ1*J1;

dA1dt = max(gammaJ1,0)*J1 + min(vA1,0)*A1 - muA1*A1;

dJ2dt = max(vA2,0)*A2 + vJ2*J2 - max(gammaJ2,0)*J2 - muJ2*J2;

dA2dt = max(gammaJ2,0)*J2 + min(vA2,0)*A2 - muA2*A2;

dYdt = [dRSdt; dRLdt; dJ1dt; dA1dt; dJ2dt; dA2dt];

end