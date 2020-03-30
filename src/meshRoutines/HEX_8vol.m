function vol = HEX_8vol(ex, ey, ez)
%HEX_8VOL Summary of this function goes here
%   Detailed explanation goes here

% Split the hexahedron into 6 tetrahedrons and sum their volumes
vol = 0;
k = 7;
i = 4;
for j = [8 3]
    points = [1, i, j, k];
    vol = vol + TET_4vol(ex(points, :), ey(points, :), ez(points, :));
end
i = 2;
for j = [6 3]
    points = [1, i, j, k];
    vol = vol + TET_4vol(ex(points, :), ey(points, :), ez(points, :));
end
i = 5;
for j = [8 6]
    points = [1, i, j, k];
    vol = vol + TET_4vol(ex(points, :), ey(points, :), ez(points, :));
end
end

