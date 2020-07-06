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
%    isFStable.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: computes spectral radius and  and returns highest eigenvalue of Jacobian of the operator g,
%    allowing the determination of the stability of a given self-orbit
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


function valP = isFStable(Me, R, V, params, DISPLAYFIGURE)

% R, V are parameters of an orbit corresponding to an equilibrum at memory Me, for
% the given experimental parameters and choice of attenuation cutOff params (computed with getExperimentalPosition())

% isFStable :  returns the largest eigenvalue of the jacobian matrix (once the
% trivial eigenvalues 1 and exp(+/- i*theta) are removed
% And assuming only 1 pair of eigenvalue might be 

% If DISPLAYFIGURE = 1 : shows plot of eigenvalues

SAVEFIGURE = 0;

%%% Experimental parameters
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

% coherence length of walker for which the field is attenuated a less than
% 1 percent
p = computeCutOff(Me);


R = R*lambda;
V = V*lambda/Tf;

% omega = V/R;
% theta = omega*Tf;
theta = (2*(asin(V*Tf/(2*R))));
% omega = (theta/Tf);

% 2x(2) rotation matrix acting on one component of state vector angle -theta
r2 = (([[cos(-theta), -sin(-theta)] ; [sin(-theta), cos(-theta)]]));
% Block diagonal rotation matrix
rot2 = kron(eye(p), r2);

x=zeros(p,1);
for ii=1:p
    x(ii) = (R*cos((1-ii)*theta));
end

y=zeros(p,1);
for ii=1:p
    y(ii) =(R*sin((1-ii)*theta));
end

jac = jacob(x, y, Me, params );

res = rot2*jac;
vp = eig(res);

% Dirty trick : works only because we assume only pair of eigenvalue might
% be greater than one, which is true in the range studied
vpp = sort (vp);
if round(abs(vpp(end)), 5) > 1
    spectralRadius = round(abs(vpp(end)), 5);
    valP = vpp(end);
else
    vpp = vpp(1:end-3);
    spectralRadius = round(abs(vpp(end)), 5);
    valP = vpp(end);
end


if (DISPLAYFIGURE)
    % zoom faible
    figure;
    hold on;
    circle(0,0,1, 1,'.r')
    plot(vp, 'k.', 'markerSize', 4);
    grid on;
    axis ([-1, 1, -1, 1]);
    rectangle('Position',[0.85 -0.13 0.27 0.26])
    
    axis equal
    box on;
    hold off
    
    
    titre = strcat(['Me=', num2str(Me), ' R=', num2str(R/lambda), ' V=' , num2str(V*Tf/lambda ), ' Spectral radius=', num2str(spectralRadius, 10)]);
    title(titre, 'interpreter', 'latex');
    
    ax = axes('Position',[0.35 0.35 0.35 0.35]);
    
    subplot(ax)
    hold on;
    circle(0,0,1, 1,'.r')
    plot(vp, 'k.');
    grid on;
    axis ([0.85, 1.12, -0.13, 0.13]);
    axis manual;
    plot(1, 0, 'ob', 'markerSize', 7, 'linewidth', 1.5);
    plot(exp(1i*theta), 'oc', 'markerSize', 7, 'linewidth', 1.5);
    plot(exp(-1i*theta), 'oc', 'markerSize', 7, 'linewidth', 1.5);
    box on;
    hold off
    % gcf.WindowState = 'maximized';
    drawnow;
    
    
    %
    if(SAVEFIGURE)
    figFileName = titre;
    saveas(gcf, [figFileName, '.fig']);
    % Export a figure in HD
    width = 120; % cm
    height = 80; % cm
    %set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54]);
    set(gcf, 'PaperPositionMode', 'auto');
    print ( [figFileName, '.tiff'], '-dtiff', '-r400');
    end
end

end

