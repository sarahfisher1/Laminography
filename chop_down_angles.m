% CHOP_DOWN_ANGLES
% -------------------------------------------------------------------------
% Function used to reduce the number of projection angles by sampling.
% Use if you need to reduce the size of your data set for testing or
% performance reasons. Use in conjunction with chop_down_sinogram to reduce
% the number of projections in the corresponding sinogram.
% 
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a 
% 
% INPUT:
%   angles - An array containing the projection angles (in order)
%   factor - The sampling factor
% 
% OUTPUT:
%   angles_chopped - The sampled angles array
% 
% EXAMPLE:
%   If the input angles array is [0, 1, 2, 3, 4, 5] and the sampling factor
%   is 2 then the output array is [1, 3, 5].
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
% -------------------------------------------------------------------------
%
% Any questions and comments please contact sarah.fisher-1@outlook.com
% -------------------------------------------------------------------------


function angles_chopped = chop_down_angles(angles, factor)

proj_num = length(angles);
k = 0;
angles_chopped = zeros(floor(proj_num/factor), 1);

for i = 1:proj_num
    if ( i/factor == double(int32(i/factor)) )
        k = k + 1;
        angles_chopped(k) = angles(i);
    end
end

end

