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
%    newPosition.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: This script is used to compute and plot the stability diagram
%    For a given experimental situatin
%       1 - computes radius and speed of possible self orbits at a given memory
%       2 - Check if this orbit is stable with linear stability analysis
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


clear all;
close all;
clc;

SAVEFIGURE = 0;
params = getParameters(1, 0.001); % This is the experimental situation described in the paper "stable self orbit ..."

f = params(1);
Tf = params(2);
m = params(3);
lambda = params(4);
fps = params(5);
kf = params(6);
Ch0 = params(7);
D = params(8);
I = params(9);
J = params(10);
K = params(11);

MeList = 40:1:145; % 49 - 200
resolution = 120;
tic;
pos = [];
for ii = 1:length(MeList)
MeList(ii)
pos =  [pos, getEquilibriumPosition(MeList(ii), 0.1, 5.1, 0.01, 0.02, resolution, params,0 ) ];
toc
end


% verif
% figure ; scatter(pos(1,:), pos(2,:));

% Examples of stability tests
% MM = pos(1,43);
% RR = pos(2,43);
% VV = pos(3,43);
% isFStable(MM, RR, VV);
% % MM = pos(1,2);
% % RR = pos(2,2);
% % VV = pos(3,2);
% % isFStable(MM, RR, VV);

% Add stability criterion to data
for iii = 1:length(pos(1,:))
    iii/length(pos(1,:))
    if (~isempty(pos(2,iii)))
    MM = pos(1,iii);
    RR = pos(2,iii);
    VV = pos(3,iii);
    tmp = isFStable(MM, RR, VV, params);
    pos(4,iii) = real(tmp);
    pos(5,iii) = imag(tmp);
    pos (6, iii ) = norm(tmp);
     end
end

% Final plot (horizontal) 
figure; hold on
colormap parula
pointsize = 30;
maxPos = length(pos(1,:)) ;
maxMem = pos( 1, maxPos);
for iii = 1:maxPos; %length(pos(1,:)) % 470 = ME145
    tmp = pos (6, iii );
    if pos(6,iii) >1
        if pos(5,iii) > 0
            scatter(pos(1,iii), pos(2,iii),  pointsize,  tmp, 'filled');
%             scatter(pos(2,iii), pos(1,iii),  40, 'or');
        else
             scatter(pos(1,iii), pos(2,iii),   pointsize,  tmp, 'filled');
%             scatter(pos(2,iii), pos(1,iii),  40, 'xr');
        end
    else
         scatter(pos(1,iii), pos(2,iii),   pointsize/4,  'k', 'filled');
%         scatter(pos(2,iii), pos(1,iii), 110, 'ro');
    end
    drawnow
end
cb = colorbar;
cb.Label.String = 'Most unstable eigenvalue modulus';
cb.Label.FontSize= 20;

% line([0 0.7], [70,70] , 'color', 'k', 'LineStyle', '--')
% line([0 1.2], [90,90] , 'color', 'k', 'LineStyle', '--')
% ax = gca;
% ax.YTick = 40:10:200;
kf = 2*pi/lambda;
grid off;
ylabel('$r/\lambda$', 'interpreter', 'latex', 'FontSize', 20);
xlabel('Memory', 'interpreter', 'latex', 'FontSize', 20);
title('Stability diagram', 'interpreter', 'latex', 'FontSize', 20);
z0 = besselzero(0,9)/(kf*lambda);
z1 = besselzero(1,9)/(kf*lambda);
for ii = 1:9
line([20, maxMem], [z0(ii), z0(ii)] , 'color', 'r', 'LineStyle', '--');
end
for ii = 1:8
line([20, maxMem], [z1(ii), z1(ii)] , 'color', 'k', 'LineStyle', '--');
end

% find(pos(1,:) == 140)
% pos(:,find(pos(1,:) == 145))
plot (140, pos(2,1437),  'xr', 'markersize', 20, 'linewidth', 1.5)






% iii = 3
% 
%     MM = pos(1,iii);
%     RR = pos(2,iii);
%     VV = pos(3,iii);
%     tmp = isFStable(MM, RR, VV);
%     pos(4,iii) = real(tmp);
%     pos(5,iii) = imag(tmp);
%     pos (6, iii ) = norm(tmp);
