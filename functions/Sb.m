function [ slopebreak ] = Sb(Z, cellsize, near_distance, far_distance, direction )
%UNTITLED Adam Winstral's slopebreak term, Sb, where breaks from the
%prevailing wind are defined as:

% the near field Sx minus the far field Sx
% 
%   
% 
%  
% 
% 
%     
% 
Sx_near = Sx(Z, cellsize, near_distance(1), near_distance(2), direction);
Sx_far = Sx(Z, cellsize, far_distance(1), far_distance(2), direction);
slopebreak = Sx_near - Sx_far;

end

