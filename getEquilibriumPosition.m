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
%    getEquilibriumPosition.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: Computes parameter of circular auto-orbit given Memory and
%    experimental parameters
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */

function res = getEquilibriumPosition(Me, rmin, rmax, vmin, vmax, resolution, params, DISPLAYFIGURE) 
% find self orbit equilibrium position in (r,v) plane 
% Search range [rmin, rmax] , [vmin, vmax] and resolution parameter
% "resolution"
% The model parameters are in params (run params = getParameters() first;
% DISPLAYFIGURE = 1 displays the figure with curves intersections

SAVEFIGURE = 0;
MAKEMOVIE = 0;

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

MeMin =Me;
MeMax = Me;

memList = MeMin:20:MeMax;

% For reference only
% r0 = 0.408*lambda;
% v0 = 0.0144*lambda/Tf;
% theta0 = 2*asin(v0*Tf/(2*r0));
% omega0 = theta0/Tf;
% v20 = r0*sin(theta0)/Tf;

% for graph up to ME 140
% Goutte sam 0.1 5.1 
% Goutte Ferrofluid 0. 5.5 
rSamples =resolution*3;
rmin = rmin*lambda;
rmax =rmax*lambda;
resrr = (rmax-rmin)/rSamples;
rr = rmin: resrr : rmax;


% Goutte sam 0.01 0.02
% Goutte Ferrofluid 0.03 0.04 %  pour ME 140
vSamples = resolution/3;
vmin = vmin*lambda/Tf;
vmax = vmax*lambda/Tf;
resvv = (vmax-vmin)/vSamples;
vv = vmin: resvv : vmax;

Gr = zeros(length(memList),length(rr), length(vv));
Gtheta = zeros(length(memList),length(rr), length(vv));

for ii = 1:length(memList)
    memList(ii);
    Me = memList(ii);
    cutOff = computeCutOff(Me);
    NCircle =cutOff+1;

  
    for jj = 1:length(rr)
        
%         jj/length(rr) % Compteur !
        r = rr(jj);
        
        for kk = 1:length(vv)
          
             v = vv(kk);
           
            theta = 2*asin(v*Tf/(2*r));
            omega = theta/Tf;
            ang = -omega.*(1:NCircle)/fps;
            xCircle =  r.*cos(ang);
            yCircle =  r.*sin(ang);
            
            %% INITIAL CONDITIONS
            rxini = xCircle(1:NCircle);
            ryini = yCircle(1:NCircle);
            
            rhx =  rxini(1) ;
            rhy =  ryini(1) ;
            
            for frameNumber = 1:NCircle-1
      
                rhx = [rxini(frameNumber+1), rhx];
                rhy = [ryini(frameNumber+1), rhy];

            end
            [Gx, Gy] = computeGradient(rhx(2:end), rhy(2:end), cutOff, lambda, Me);
            G = [Gx, Gy];
            ur = ([rhx(2)/lambda, rhy(2)/lambda]);
            ur = ur/norm(ur);
            Gr(ii, jj, kk) = dot(G, ur);
            
            
    %  tangential component
va = [rhx(1) - rhx(2), rhy(1) - rhy(2)]'/Tf;
vb = [rhx(2) - rhx(3), rhy(2) - rhy(3)]'/Tf;
v2 = (va+vb)/2;
uv = v2/norm(v2);
Gtheta(ii, jj, kk) = dot(G, uv);
         
        end
    end
    
end


[RR, VV] = meshgrid(rr, vv);

% figure;
for ii = 1:length(memList)

% Radial force
FwR = squeeze(-Ch0*Gr(ii,:,:));
FwR = FwR';
Fcentrifugal = inline('m.*v.^2./r', 'm', 'r','v');
Fc = Fcentrifugal(m, RR, VV);
FTotR = Fc +FwR;

% Tangent force
FwT = squeeze(-Ch0*Gtheta(ii,:,:));
FwT = FwT';

Fviscous = inline('-D.*v + 0.*r', 'D', 'r', 'v');
V2 = ( RR./Tf.*sin( 2.*asin( VV.*Tf./(2.*RR) ) ) );
Fv = Fviscous(D, RR, V2);
FTotT =  Fv + FwT;
% Theoretical value at infinite memory
% FwThetaTh = C/(v0*Tf)*(1- (besselj(0,2.*pi/lambda.*rr)).^2 );
% FTotTTh = Fv + FwThetaTh;

fig = figure ('visible','off');hold on;
[mt, ct] = contour(RR./lambda, VV.*Tf./lambda, FTotT, [0 0], 'color', 'k');
[mr, cr] = contour(RR./lambda, VV.*Tf./lambda, FTotR, [0 0], 'color', 'b');
close(fig)

ind = find (mt(2,:)>vmax*Tf/lambda | mt(2,:)< vmin*Tf/lambda);
mt(:,ind) = NaN;
ind = find (mt(1,:)>rmax/lambda | mt(1,:)< rmin/lambda);
mt(:,ind) = NaN;

ind = find (mr(2,:)>vmax*Tf/lambda | mr(2,:)< vmin*Tf/lambda);
mr(:,ind) = NaN;
ind = find (mr(1,:)>rmax/lambda | mr(1,:)< rmin/lambda);
mr(:,ind) = NaN;

% 
[x0,y0,iout,jout] = intersections(mt(1,:),mt(2,:),mr(1,:),mr(2,:), 1);
[ind] = find( y0<vmax*Tf/lambda & y0>vmin*Tf/lambda);
rrr = x0(ind);
vvv = y0(ind);


res(3,:) = vvv';
res(2,:) = rrr';
res(4,:) = 1;
res(1,:) = Me;

% 
% 
% % Radial
% figure;
% %subplot(length(memList),3,sub2ind([ 3 length(memList)], 1,ii));
% subplot(1,2,1);
% h = pcolor(RR./lambda, VV.*Tf./lambda, FTotR);
% set(h, 'EdgeColor', 'none');
% shading interp;
% axis('xy');axis square;
% title(['Total Radial Force'], 'interpreter', 'latex');
% xlabel('$r/\lambda$', 'interpreter', 'latex');
% ylabel('$v/v_{\varphi}$', 'interpreter', 'latex')
% hold on;
% [mr, cr] = contour(RR./lambda, VV.*Tf./lambda, FTotR, [0 0], 'color', 'k');
% c = colorbar;
% c.Label.String = 'Force (N)';
% %suptitle(['  cutOff = ', num2str(cutOff), '  Me = ', num2str(Me)]);
% 
% 
% % Tangent
% %subplot(length(memList),3,sub2ind([ 3 length(memList)], 2,ii));
% subplot(1,2,2);
% h = pcolor(RR./lambda, VV.*Tf./lambda, FTotT);
% set(h, 'EdgeColor', 'none');
% shading interp;
% axis('xy');axis square;
% title(['Total Tangent Force'], 'interpreter', 'latex');
% xlabel('$r/\lambda$', 'interpreter', 'latex');
% ylabel('$v/v_{\varphi}$', 'interpreter', 'latex')
% hold on;
% [mt, ct] = contour(RR./lambda, VV.*Tf./lambda, FTotT, [0 0], 'color', 'k');
% c = colorbar;
% c.Label.String = 'Force (N)';
% hold off;
% suptitle(['  cutOff = ', num2str(cutOff), '  Me = ', num2str(Me)]);
% 
% if(SAVEFIGURE)
% figFileName = strcat(['Me', num2str(Me), 'cutOff = ', num2str(cutOff), 'totalForces']);
%  saveas(gcf, [figFileName, '.fig']);
% % Export a figure in HD
% width = 120; % cm 
% height = 80; % cm
% %set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54]);
% set(gcf, 'PaperPositionMode', 'auto');
% print ( [figFileName, '.tiff'], '-dtiff', '-r400');
% end
% 
 % Intersection
 %subplot(length(memList),3,sub2ind([3, length(memList)], 3,ii));
if DISPLAYFIGURE
 figure;
hold on;
[mt, ct] = contour(RR./lambda, VV.*Tf./lambda, FTotT, [0 0], 'color', 'k');
[mr, cr] = contour(RR./lambda, VV.*Tf./lambda, FTotR, [0 0], 'color', 'k');
 plot(rrr, vvv, 'rx', 'markersize', 10, 'linewidth', 1.5);axis square;
title(['Equilibrium positions'], 'interpreter', 'latex', 'FontSize', 15);
xlabel('$r/\lambda$', 'interpreter', 'latex','FontSize', 15);
ylabel('$v/v_{\varphi}$', 'interpreter', 'latex','FontSize', 15)
axis square;
grid on;
%suptitle(['  NCircle = ', num2str(NCircle), '  Me = ', num2str(Me)]);
% 
if(SAVEFIGURE)
figFileName = strcat(['Me', num2str(Me), 'cutOff = ', num2str(cutOff), 'Intersection']);
saveas(gcf, [figFileName, '.fig']);
% Export a figure in HD
width = 120; % cm 
height = 80; % cm
%set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54]);
set(gcf, 'PaperPositionMode', 'auto');
print ( [figFileName, '.tiff'], '-dtiff', '-r400');
end
% 
% end
end

end


% verif
%figure ; scatter(pos(2,:), pos(1,:));










