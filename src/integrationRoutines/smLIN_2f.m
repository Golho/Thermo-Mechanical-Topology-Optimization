function [fe] = smLIN_2f(ex, ey, loadMagnitude, thickness)
%SMLIN_2F Summary of this function goes here
%   Detailed explanation goes here
fe = zeros(4, 1);
fe([1 3]) = hLIN_2f(ex, ey, loadMagnitude(1), thickness);
fe([2 4]) = hLIN_2f(ex, ey, loadMagnitude(2), thickness);
end

