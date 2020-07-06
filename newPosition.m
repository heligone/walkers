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
%    Purpose: Computes new position of the walker's chain of impact given its previous
%    coordinates z
%    @author Samuel BERNARD-BERNARDET
%    @version 1.0 06/07/2020
% */

function znew = newPosition(z, p, I, J, K, kf, Me)
% Computes new position of the walker's chain of impact given its previous
% coordinates z
% p number of impacts taken into account
% I, J, K are th eiterative model parameter
% kf is the faraday wave number
% Me is the Memory of the bath

znew(1) = I*z(1) + J*z(2);
for k=2:p
znew(1) = znew(1) + K*exp(-(k-1)/Me)*besselj(1, kf*norm( z(1)-z(k) ) ) / norm( z(1)-z(k) ) * (z(1)-z(k));
end

for ii = 2:length(z)
znew(ii) = z(ii-1) ;
end

end