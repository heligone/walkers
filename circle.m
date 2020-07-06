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
%    circle.m
%    Purpose: draws circle from center coordinates and radius
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */


function circle(x,y,r, size, color)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.0001 angle step
ang=0:0.0001:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp, color, 'markersize', size ,'HandleVisibility','off' );
end