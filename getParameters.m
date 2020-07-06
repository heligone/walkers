% /*
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%     http://www.apache.org/licenses/LICENSE-2.0
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
% */
% 
% /**
%    @file
%    WALKERS V1.0 project
%    computeCutOff.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: retruns experimental parameters corresponding to 3 different
%    expermiental situations
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


function params =getParameters( ExperienceNumber, attenuationAtCutoff)
%%  Returns an array with all parameters corresponding to a given experimetal situation
%%  This function must be run before any other

% ExperienceNumber = 1 --> droplet described in paper 'stable self spin of
% a wave dressed particle'

% ExperienceNumber = 2 -->  droplet described in paper 'self attraction
% into spinning eigenstates...', Fig1 (Mass is computed from the volume of
% the droplet

% ExperienceNumber = 3 -->  droplet described in paper 'self attraction
% into spinning eigenstates...' fig4, with a mass modified by ferrofluid 

%% attenuationAtcutoff specifies the level of field attenuation after which the influence of prevous impacts are neglected

% The returned parameters are as follows : 
% f = params(1); forcing freqeuncy
% Tf = params(2); Faraday frequency
% m = params(3); mass of droplet
% lambda = params(4); Faraday Wavelength
% fps = params(5); Faraday frequency
% kf = params(6); Faraday wave number
% Ch0 = params(7); Wave field coupling parapeter
% D = params(8); Viscous disspation parameter
% I = params(9); discrete iterative model parameter
% J = params(10); discret iterative model parameter
% K = params(11); discret iterative model parameter

global fieldAttenuationAtCutoff ;
 fieldAttenuationAtCutoff = attenuationAtCutoff;

%%% Memory for analysed circular experimental orbit

if ExperienceNumber == 1
%%%%%%%%%%%%%%%%%%%%%%%%%
% SBB Droplet
Me_exp = 70;
%%% Droplet diameter m
diameter = 460*10^-6; 
%%% bath shaking freq
f = 70 ; 
Tf = 2/f;
%%% Faraday WaveLength
lambda = 5.2*10^-3;
%  Radius and speed of experimental orbit used for analysis
r_exp = 0.408*lambda;
v_exp = 0.0144*lambda/Tf;
%%%%%%%%%%%%%%%%%%%%%%%%%
end


if (ExperienceNumber == 2) | (ExperienceNumber == 3)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Emmanuel droplet
Me_exp = 140;
%%% Droplet diameter m
diameter = 700*10^-6; 
%%% bath shaking freq
f = 80 ; 
Tf = 2/f;
%%% Faraday WaveLength
lambda = 4.75*10^-3;
%  Radius and speed of experimental orbit used for analysis
r_exp = 0.385*lambda;
v_exp = 0.036842*lambda/Tf;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% Oil specifications
sigma = 0.0206;  % surface tension, N/m
rho = 950; % density, Kg/m^3
% nu_dynamic = 20*10^(-3); % oil 20 Cst
% nu = nu_dynamic/rho ; % m²/s kinematic viscosity
nu = 20*10^(-6); % oil 20 Cst m/s

% Air specifications
% from 'drops walking on a vibrating bath'
rho_a = 1.2;%  kg/m3
mu_a = 1.84*10^-5;

g = 9.81; % gravity, m/s^2

 % impact neglected if amplidtude after p*Tf < 0.001*h0; inf for no cutoff ; 
p_exp = computeCutOff(Me_exp);
%p = cutOff; % Just another name for this variable for compatibility with other scripts

% Geometry
R_Drop = diameter/2; % Droplet diameter m
volume = 4/3*pi*(R_Drop)^3 ;% m3
m = rho.*volume;  %kg

%%%%%
if (ExperienceNumber == 3)
% Override mass calculation from volum in case of ferrofluid 
%(Efective mass)
 m = 3.82e-7 ; % goutte emmanuel masse efefctive ferrofluide
%%%%%
end

% Faraday Wave
fps = f/2;% Faraday frequency = camera frequency
kf = 2*pi/lambda;

%% Coupling constant from from "hydrodynamic spin states" for information and comparison  only
C_Drag = 0.17; % drag coefficient 
sin_phi = 0.2; % Phase at impact
D_Hydro = C_Drag*m*g*sqrt(rho*R_Drop/sigma) + 6*pi*mu_a*R_Drop*(1+rho_a*g*R_Drop/(12*mu_a*f)); % kg/s
A = sqrt(8*pi*nu*Tf)/3*(kf*R_Drop)^3*sin_phi/( 3*kf^2*sigma/(rho*g)+1 ); % Amplitude of waves m 
Ch0_Hydro = m*g*A % N.m
D_Hydro

%% coupling constant from simulation ion circular trajectory

NCircle =2000; % Number of impacts on initial circle


% Polygon calculation of angle between 2 impacts
theta_exp = 2*asin(v_exp*Tf/(2*r_exp));
omega_exp = theta_exp/Tf;

%% INITIAL circular CONDITIONS
ang = -omega_exp.*(1:NCircle)/fps; % clockwise motion
xCircle =  r_exp.*cos(ang);
yCircle =  r_exp.*sin(ang);

rxini = xCircle(1:NCircle);
ryini = yCircle(1:NCircle);

rhx =  rxini(1) ;
rhy =  ryini(1) ;

for frameNumber = 1:NCircle-1
    
    rhx = [rxini(frameNumber+1), rhx];
    rhy = [ryini(frameNumber+1), rhy];
    
end

% % Constant calculations after long circle
[Gx, Gy] = computeGradient(rhx(2:end), rhy(2:end), p_exp, lambda, Me_exp);
G = [Gx, Gy];
% Radial component
ur = ([rhx(2)/lambda, rhy(2)/lambda]);
ur = ur/norm(ur);
Gr = dot(G, ur);
Gr;
GrVect = Gr.*ur;
F = m*v_exp^2/r_exp % Radial force on the droplet
CC = F/(Gr);
% alternative computation with dicret polyfon instead of circle
%  FTest = m*4*r_exp/Tf^2*(sin(theta_exp/2))^2 
% CCTest = FTest/(Gr)

% tangential component
utheta = [-rhy(2)/lambda, rhx(2)/lambda];
utheta = utheta/norm(utheta);
va = [rhx(1) - rhx(2), rhy(1) - rhy(2)]'/Tf;
vb = [rhx(2) - rhx(3), rhy(2) - rhy(3)]'/Tf;
vv = (va+vb)/2;
uv = vv/norm(vv);
Gtheta = -dot(G, uv);
Gthetavect = Gtheta.*uv;

DD = CC*Gtheta/norm(vv);

% 2 ways of calculating v2, the instant speed at impact
%(mean of speed before/after impact
% from simulation
% norm(vv)
% from polygon geometry
% v2 = r*sin(theta)/Tf;

Ch0 = CC
D = DD

% % Mathieu's value
% Ch0 = 2.32e-11
% D = 7.8e-6
%%m = 3.82e-7; % meff described in the paper

% Discrete simulation constants
I = 1 + ( 1-D*Tf/(2*m) ) / (1+D*Tf/(2*m)) ;
J = - (1-D*Tf/(2*m)) / (1+D*Tf/(2*m)) ;
K = Ch0*Tf^2*kf/(( 1+D*Tf/(2*m) ) * m); 


params = [f, Tf, m, lambda, fps, kf,  Ch0, D, I, J, K,];



end

% f = params(1);
% Tf = params(2);
% m = params(3);
% lambda = params(4);
% fps = params(5);
% kf = params(6);
% Ch0 = params(7);
% D = params(8);
% I = params(9);
% J = params(10);
% K = params(11);

