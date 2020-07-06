
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
%    jacob.m
%    Code associated to "Stable spin mode of a wave dressed particle"
%    Purpose: % Computes jacobian of operator f at point x, y (See Annex), at Memory Me
% and for experimental parameters params
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */

function jac = jacob(x,y, Me, params)
% Computes jacobian of operator f at point x, y (See Annex), at Memory Me
% and for experimental parameters params

% Experimental constant
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
% Number of impacts taken into accounts
p = computeCutOff(Me);


jac = zeros(2*p, 2*p);

jac(1,1) = I;

for k=2:p
    jac(1,1) = jac(1,1) +  (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2)))/((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*((besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*x(1) - 2*x(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)) - (kf*besselj(0, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*x(1) - 2*x(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2)))*(x(1)-x(k)))/((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*x(1) - 2*x(k))*(x(1)-x(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(3/2));

end

%jac(1,2) = 0;
for k=2:p
    jac(1,2) = jac(1,2)   - (K*exp(-(k-1)/Me)*((besselj(1, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)) - (kf*besselj(0, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2)))*(x(1) - x(k)))/((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k))*(x(1) - x(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(3/2));
end

jac(1,3) = J - (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*((besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*x(1) - 2*x(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)) - (kf*besselj(0, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*x(1) - 2*x(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))*(x(1) - x(2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*x(1) - 2*x(2))*(x(1) - x(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(3/2));

jac(1,4) = (K*exp(-1/Me)*((besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)) - (kf*besselj(0, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))*(x(1) - x(2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2))*(x(1) - x(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(3/2));

for jj= 3:p
    jac(1, 2*jj-1) = (K*exp(-(jj-1)/Me)*((besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)) - (kf*besselj(0, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))*(x(1) - x(jj)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) - (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) + (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj))*(x(1) - x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(3/2));
    jac(1, 2*jj) = (K*exp(-(jj-1)/Me)*((besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)) - (kf*besselj(0, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))*(x(1) - x(jj)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) + (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj))*(x(1) - x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(3/2));
end

%jac(2,1) = 0;
for k=2:p
    jac(2,1) = jac(2,1)   - (K*exp(-(k-1)/Me)*((besselj(1, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)) - (kf*besselj(0, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2)))*(x(1) - x(k)))/((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(1/2))*(2*y(1) - 2*y(k))*(x(1) - x(k)))/(2*((x(1) - x(k))^2 + (y(1) - y(k))^2)^(3/2));
end

jac(2,2) = I;
for k=2:p
    jac(2,2) = jac(2,2) +  (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2)))/((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*((besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)) - (kf*besselj(0, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*y(1) - 2*y(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2)))*(y(1)-y(k)))/((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2) - (K*exp(-(k-1)/Me)*besselj(1, kf*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(1/2))*(2*y(1) - 2*y(k))*(y(1)-y(k)))/(2*((x(1)-x(k))^2 + (y(1)-y(k))^2)^(3/2));
end

jac(2,4) = J - (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*((besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)) - (kf*besselj(0, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))*(y(1) - y(2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2))*(y(1) - y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(3/2));

jac(2,3) = (K*exp(-1/Me)*((besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)) - (kf*besselj(0, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2)))*(x(1) - x(2)))/((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2) + (K*exp(-1/Me)*besselj(1, kf*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(1/2))*(2*y(1) - 2*y(2))*(x(1) - x(2)))/(2*((x(1) - x(2))^2 + (y(1) - y(2))^2)^(3/2));

for jj= 3:p
    jac(2, 2*jj) = (K*exp(-(jj-1)/Me)*((besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)) - (kf*besselj(0, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))*(y(1) - y(jj)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) - (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) + (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*y(1) - 2*y(jj))*(y(1) - y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(3/2));
    jac(2, 2*jj-1) = (K*exp(-(jj-1)/Me)*((besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)) - (kf*besselj(0, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2)))*(y(1) - y(jj)))/((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2) + (K*exp(-(jj-1)/Me)*besselj(1, kf*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(1/2))*(2*x(1) - 2*x(jj))*(y(1) - y(jj)))/(2*((x(1) - x(jj))^2 + (y(1) - y(jj))^2)^(3/2));
end

for i = 3:2*p
    jac(i,i-2) = 1;
end

end


