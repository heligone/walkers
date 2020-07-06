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
%    runSimul.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: Runs simulation of walker of memory Me, with an initial circular trajectory of radius R
%    at speed V. Experimental parameters are computed first with getParameters() and stored in params variable. 
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


function runSimul(Me, R, V, params , varargin)

%   Runs simulation of walker of memory Me, with an initial circular trajectory of radius R
%   at speed V. Experimental parameters are computed first with getParameters() and stored in params variable. 
% 
%   5 optionnal inputs with defaults values
%       L=2000 Impacts in initial circular trajectory
%       NFramesimul=15000 duration of simulation
%       DISPLAY = 1 shows LIVE simulation , else, computes and plots full trajectory
%       FIELD_DISPLAY = 0 ; shows height field if FIELD_DISPLAY = 1
%       SCHLIEREN = 0 ; shows field in a schlieren style if SCHLIEREN = 1

%%% We only want 5 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 5
    error('4 optionnal inputs only');
end
optargs = {2000 15000 1 0 0}; % Set default values
optargs(1:numvarargs) = varargin;
[L, NFrameSimul, DISPLAY, FIELD_DISPLAY, SCHLIEREN] = optargs{:};
%%%

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
%%%
%%% Dimension of simulation = number of impact taken int account at 
p = computeCutOff(Me);
%%%

% Other simulation display parameters
% Display scale
scale = 3*lambda;
SET_IMAGESC = SCHLIEREN; % Simulation of schlieren light
SET_PCOLOR  = 1-SET_IMAGESC; % Height field view
%%%


R = R*lambda;
V = V*lambda/Tf;

theta = 2*asin(V*Tf/(2*R));
% omega = theta/Tf;

if (FIELD_DISPLAY)
% Cyclic colormap
hmap(1:256,1) = linspace(0,1,256);
%hmap(1:256,1) = [linspace(0.5,1,128) linspace(0,0.5,128)];
hmap(:,[2 3]) = 0.7; %brightness
huemap = hsv2rgb(hmap);
%huemap = hsv2rgb(hsv(256))
% colormap(huemap);


% siz = 2*lambda; % Size of square
siz = scale;
XMin = -siz;
XMax = siz;
YMin = -siz;
YMax = siz;
NPix = 80; %Number of pixels along dimension
if SCHLIEREN
    NPix = 200;
end

[X, Y] = meshgrid(XMin:(XMax-XMin)/(NPix-1):XMax, YMin:(XMax-XMin)/(NPix-1):YMax);
X = X./lambda;
Y = Y./lambda;
Z = zeros(size(X));

end


% Initial conditions : on a circle of radius R at speed V
z0 = zeros(p, 1);
for ii = 1:p
    z0(ii) = R*exp(-(ii-1)*1i*theta);
    if (FIELD_DISPLAY)
        Z=Z.*exp(-1/Me) + besselj(0, 2*pi/lambda*sqrt( (lambda.*X-real(z0(p-ii+1)) ).^2 + ( lambda.*Y-imag(z0(p-ii+1))).^2 ) );
    end
end
z=zeros(p, NFrameSimul);
z(:,1) = z0(1:p);

% scale = 1.5/lambda;

if DISPLAY
    figure;    
    if (FIELD_DISPLAY)
         if (SET_IMAGESC); colormap(huemap); axis xy; end
         if (SET_PCOLOR); colormap(jet);end
        axis([-scale/lambda, scale/lambda, -scale/lambda, scale/lambda]);
        axis square;
        
        if (SET_IMAGESC)
            [FX,FY] = gradient(Z);
            [theta, pho] = cart2pol (FX, FY);
            
            h = imagesc( -theta+pi, 'alphadata', 2*abs(pho)/max(abs(pho(:))));
        h = imagesc( [XMin/lambda, XMax/lambda], [XMin/lambda, XMax/lambda],-theta+pi, 'alphadata', 0.9*pho);
        
        axis xy;
        end
        
        if (SET_PCOLOR)
        h = pcolor( X, Y, Z);
        set(h, 'EdgeColor', 'none');
        shading interp; 
        end
        
    end
    
     if (FIELD_DISPLAY); hold on;end
    plot(z(1:min([L,p,Me]),1)/lambda, '.', 'markersize', 6);
    axis([-scale/lambda, scale/lambda, -scale/lambda, scale/lambda]);
    axis square;
    grid on;
      if (FIELD_DISPLAY) hold off;end 
    title(ii);
    drawnow;

end

% Case of initial conditions shorter than p :
% Only impacts during initial conditions  taken into account
for ii=2:p-L+1
    z(:,ii) = newPosition(z(:, ii-1), L+ii-2, I, J, K, kf, Me);
    if (FIELD_DISPLAY)
          Z=Z.*exp(-1/Me) + besselj(0, 2*pi/lambda*sqrt( (lambda.*X-real(z(1,ii)) ).^2 + ( lambda.*Y-imag(z(1,ii))).^2 ) );
    end

    if DISPLAY
        if (FIELD_DISPLAY)
        if (SET_IMAGESC)
            [FX,FY] = gradient(Z);
            [theta, pho] = cart2pol (FX, FY);
           
            h = imagesc( -theta+pi, 'alphadata', 2*abs(pho)/max(abs(pho(:))));
        h = imagesc( [XMin/lambda, XMax/lambda], [XMin/lambda, XMax/lambda],-theta+pi, 'alphadata', 0.9*pho);
       axis xy;
        end
        
        if (SET_PCOLOR)
        h = pcolor( X, Y, Z);
        set(h, 'EdgeColor', 'none');
        shading interp; 
        end
            
        end
         if (FIELD_DISPLAY); hold on;end
        plot(z(1:Me,ii)/lambda, '.', 'markersize', 6);
            axis([-scale/lambda, scale/lambda, -scale/lambda, scale/lambda]);
        axis square;
        grid on;
          if (FIELD_DISPLAY); hold off;end
        title(ii);
        drawnow;
       
    end
end

% Rest of trajectory after initial conditions has been processed
for ii= max([p-L+1 2]):NFrameSimul
    ii/NFrameSimul
    z(:,ii) = newPosition(z(:, ii-1), p, I, J, K, kf, Me);
        if (FIELD_DISPLAY)
          Z=Z.*exp(-1/Me) + besselj(0, 2*pi/lambda*sqrt( (lambda.*X-real(z(1,ii)) ).^2 + ( lambda.*Y-imag(z(1,ii))).^2 ) );
        end
    
    if DISPLAY
            if (FIELD_DISPLAY)
        if (SET_IMAGESC)
            [FX,FY] = gradient(Z);
            [theta, pho] = cart2pol (FX, FY);
         
            h = imagesc( -theta+pi, 'alphadata', 2*abs(pho)/max(abs(pho(:))));
        h = imagesc( [XMin/lambda, XMax/lambda], [XMin/lambda, XMax/lambda],-theta+pi, 'alphadata', 0.9*pho);
        axis xy;
        end
        
        if (SET_PCOLOR)
        h = pcolor( X, Y, Z);
        set(h, 'EdgeColor', 'none');
        shading interp; 
        end
       
            end
         if (FIELD_DISPLAY); hold on;end
        plot(z(1:Me,ii)/lambda, '.', 'markersize', 6);
            axis([-scale/lambda, scale/lambda, -scale/lambda, scale/lambda]);
        axis square;
        grid on;
        title(ii);
         if (FIELD_DISPLAY); hold off;end
        drawnow;
        
    end
    
    
end

% Full trajectory
figure;
plot(z(1,1:NFrameSimul)/lambda);
axis equal;
grid on;

end