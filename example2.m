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
%    example2.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: script showing workflow for analysing Experimental data corresponding to auto-orbit described in "self-attraction into spinning eigenstates.."
%    The mass of the droplet is the effefctive mass due to the presence of  ferrofluid in the droplet (fig 4)
%   with Cutoff at an attenuation of 0.1 %
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


%%%%%%%%%%%
% Run FIRST
%%%%%%%%%%%

close all;
clear all;

% Experimental data corresponding to auto-orbit described in "self-attraction into spinning eigenstates.."
% The mass of the droplet is the effefctive mass due to the presence of  ferrofluid in the droplet (fig 4)
% with Cutoff at an attenuation of 0.1 %
params = getParameters(3, 0.001);

%%%%%%%%%%%%%%%
% END Run FIRST
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each of the following line can be run independantly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Me = 140 Unstable and no bound state
%%%%%%%%%
getEquilibriumPosition(140, 0.1, 5.1, 0.03, 0.04,  120, params, 1) % Zoom out
getEquilibriumPosition(140, 0.3, 0.4, 0.033, 0.038,  120, params, 1) % Zoom in

% Point A
isFStable(140,  0.385, 0.036842, params, 1)
% Simulation showing droplet trajectory and its coherence tail
runSimul(140, 0.385, 0.036842, params, 2000, 15000, 1, 0, 0) % Instability appears at iteration 3500
% Simulation showing droplet trajectory and its coherence tail, plus the  height field
runSimul(140, 0.385, 0.036842, params, 2000, 15000, 1,1, 0 )
% Simulation showing droplet trajectory and its coherence tail, plus the  gradient of the wave field
% The color endodes the direction of the gradient, the transparency indicates the norm 
% (White means a small gradient)
runSimul(140, 0.385, 0.036842, params, 2000, 15000, 1,1, 1 )

%%%%%%%%%
% Me = 50 Unstable but bound state
%%%%%%%%%
getEquilibriumPosition(50, 0.1, 5.1, 0.03, 0.04,  120, params, 1) % Zoom out
getEquilibriumPosition(140, 0.3, 0.4, 0.033, 0.038,  120, params, 1) % Zoom in

% Point A
isFStable(50, 0.392766, 0.036350, params, 1)
% Simulation showing droplet trajectory and its coherence tail
runSimul(50, 0.392766, 0.036350, params, 2000, 15000, 1, 0, 0) % Instability appears at iteration 3500
% Simulation showing droplet trajectory and its coherence tail, plus the height field
runSimul(50, 0.392766, 0.036350, params, 2000, 15000, 1,1, 0 )

%%%%%%%%%
% Me = 25 stable
%%%%%%%%%
getEquilibriumPosition(25, 0.1, 5.1, 0.03, 0.04,  120, params, 1) % Zoom out
getEquilibriumPosition(25, 0.3, 0.4, 0.033, 0.038,  120, params, 1) % Zoom in

% Point A
isFStable(25, 0.42344, 0.034549, params, 1)
% Simulation showing droplet trajectory and its coherence tail
runSimul(25, 0.42344, 0.034549, params, 2000, 15000, 1, 0, 0) 
% Simulation showing droplet trajectory and its coherence tail, plus the height field
runSimul(25, 0.42344, 0.034549, params, 2000, 15000, 1, 1, 0 )


