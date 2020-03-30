function vol = TET_4vol(ex, ey, ez)
%TET_4VOL Summary of this function goes here
%   Detailed explanation goes here
% Translate the element so that point 1 coincides with origin
ex(:, :) = ex(:, :) - ex(1, :);
ey(:, :) = ey(:, :) - ey(1, :);
ez(:, :) = ez(:, :) - ez(1, :);
vols = ex([4 3 4 2 3 2], :).*ey([3 4 2 4 2 3], :).*ez([2 2 3 3 4 4], :);
vol = 1/6*abs(-vols(1, :)+vols(2, :)+vols(3, :)-vols(4, :)-vols(5, :)+vols(6, :));
end

