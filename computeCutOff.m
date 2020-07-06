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
%    Purpose: Computes number of impacts after which their influence is to
%    be neglected at a given Memory and at a given chosen attenuation
%    thresold  (global fieldAttenuationAtCutoff)
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


function p = computeCutOff(Me)
% At a given memory Me, we neglect the influence of impacts that occured p
% Faraday periods before the current impact.

% p is determined by the following condition : the influence of the impact that
% occured p Faraday period before current impact is
% fieldAttenuationAtCutoff times  the maximum field created by
% a single impact
% For example : fieldAttenuationAtCutoff = 0.001
% The global variable fieldAttenuationAtCutoff comes form getPArameters
% function that must be called first.


global fieldAttenuationAtCutoff ;
p = floor(-Me*log(fieldAttenuationAtCutoff));


end

