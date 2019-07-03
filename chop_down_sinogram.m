% CHOP_DOWN_SINOGRAM
% -------------------------------------------------------------------------
% Function used to reduce the number of projection images in a sinogram by
% sampling. Use if you need to reduce the size of your data set for testing
% or performance reasons. Use in conjuction with chop_down_angles to reduce
% the number of projection angles in the corresponding projection angles list. 
% 
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a
% 
% INPUT:
%   sinogram - a matrix containing the sinogram you want to sample from
%   factor - the sampling factor 
% 
% OUTPUT:
%   sinogram_chopped - the sampled sinogram
%
% -------------------------------------------------------------------------
% Copyright (c) 2019 S L Fisher D J Holmes J S Jørgensen P Gajjar J Behnsen W R B Lionheart P J Withers
% University of Manchester
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Contact: sarah.fisher-1@outlook.com
% -------------------------------------------------------------------------

function sinogram_chopped = chop_down_sinogram_v2(sinogram, factor)

% Take a sample of every 10 angles 

N_1 = length(sinogram(:,1,1));
N_2 = length(sinogram(1,1,:));
proj_num = length(sinogram(1,:,1));
k = 0;


sinogram_chopped = zeros(N_1, floor(proj_num/factor), N_2);

for i = 1:proj_num
   
    if ( i/factor == double(int32(i/factor)) )
        k = k + 1;
        sinogram_chopped(:, k, :) = sinogram(:, i, :);
        
    end
    
end


end

