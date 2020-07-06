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
% /*
%    @file
%    WALKERS V1.0 project
%    example1.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: script showing workflow for analysing Experimental data corresponding to auto-orbit described in  auto-orbit described in "stable spin mode of a wave-dressed particle"
%    with Cutoff at an attenuation of 0.1 %
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020


%%%%%%%%%%%
% Run FIRST
%%%%%%%%%%%

close all;
clear all;

% Experimental data corresponding to situation 1 : auto-orbit described in 'stable spin mode of a wave-dressed particle"
% with Cutoff at an attenuation of 0.001 = 0.1 % (That parameter can be played with)
params = getParameters(1, 0.001);

%%%%%%%%%%%%%%%
% END Run FIRST
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each of the following line can be run independantly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Me = 70 
%%%%%%%%%
% Graphical determination of possible orbits
getEquilibriumPosition(70, 0.1, 5.1, 0.01, 0.02,  120, params, 1); % Zoom out
getEquilibriumPosition(70, 0.4, 0.6, 0.01, 0.02,  120, params, 1) % Zoom in

% The following simulatiosn use orbit parameters computed with Cutoff at an attenuation of 0.001
% If the attenuation at cutOff is changed, the equilibrium positions must  be updated with getEquilibriumPosition() 

% Point A
% Is this orbit stable (Linear stability eigenvalue analysis)
isFStable(70,  0.408, 0.0144, params, 1)

% Simulation showing droplet trajectory and its coherence tail ( 2000 impacts on initial circular trajectory)
runSimul(70, 0.408, 0.0144, params, 2000, 15000, 1, 0, 0) %
% Simulation showing droplet trajectory and its coherence tail, plus the height field
runSimul(70, 0.408, 0.0144, params, 2000, 15000, 1,1, 0 )
% Simulation showing droplet trajectory and its coherence tail, plus the gradient of the wave field
% the color endodes the direction of the gradient, the transparency indicates the norm of the gradient.
% (White means a small gradient)
runSimul(70, 0.408, 0.0144, params, 2000, 15000, 1,1, 1 )

% point B

isFStable(70, 0.581132729, 0.013419207, params, 1)
runSimul(70, 0.581132729, 0.013419207, params, 2000, 15000, 1, 1)

%%%%%%%%%
% Me = 90
%%%%%%%%%
getEquilibriumPosition(90, 0.4, 0.6, 0.01, 0.02,  120, params, 1)

% Point C Stable
isFStable(90, 0.39757839, 0.014684946, params, 1)
runSimul(90, 0.39757839, 0.014684946, params,  2000, 15000, 1, 1, 0)

% Point D radially unstable
isFStable(90, 0.594444, 0.013568, params, 1)
runSimul(90, 0.594444, 0.013568, params,  2000, 15000, 1, 1)

% Point E oscillatory unstable
isFStable(90, 0.921109, 0.014301, params, 1);
runSimul(90, 0.921109, 0.014301, params,  2000, 15000, 1, 1, 0)
% Plots whole trajectory, showing instability and attraction into a smaller orbit
runSimul(90, 0.921109, 0.014301, params,  2000, 5000, 0, 0, 0) 

% Point F Radially unstable
isFStable(90, 1.07317, 0.014008, params, 1);
runSimul(90, 1.07317, 0.014008, params,  2000, 15000, 1, 1, 0)
runSimul(90, 1.07317, 0.014008, params,  2000, 5000, 0, 0, 0) % Plot final trajectory
