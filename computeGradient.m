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
%    computeGradient.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: computes field gradient at given coordinate, given history of
%    bounces, threshold, wavelength and Memory
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */

function [Gx, Gy] = computeGradient(rhx, rhy, p, lambda, Me)
%   computes gradient of field at rhx(1), rhy(1) considering previous impacts stored in rhx and rhy (From indices 2 onwards)  
%   p is a cutoff : the maximum number of impacts taken into account

rx = rhx(1);
ry = rhy(1);
Gx = 0;
Gy = 0;
endLoop = min(p, length(rhx));
for ii = 2:endLoop
    impact = ii;
    Gx = Gx - exp(-(ii-1)/Me)* (2*pi/lambda).*besselj(1,  (2*pi/lambda).*sqrt( (rx-rhx(impact)).^2 + (ry-rhy(impact)).^2 ) ) .* (rx-rhx(impact)) ./ sqrt( (rx-rhx(impact)).^2 + (ry-rhy(impact)).^2);
    Gy = Gy - exp(-(ii-1)/Me)* (2*pi/lambda).*besselj(1,  (2*pi/lambda).*sqrt( (rx-rhx(impact)).^2 + (ry-rhy(impact)).^2 ) ) .* (ry-rhy(impact)) ./ sqrt( (rx-rhx(impact)).^2 + (ry-rhy(impact)).^2);
end
 
end

