function [m] = rpy2matrix(xyz)
% XYZ rotation
% Function to return a generalised rotation matrix (4x4) given an input set
% of Euler angles.

%Written by Matthew Field

c1 = cos(xyz(1)); c2 = cos(xyz(2)); c3 = cos(xyz(3));
s1 = sin(xyz(1)); s2 = sin(xyz(2)); s3 = sin(xyz(3));

[m] = [[c3*c2   c3*s2*s1-s3*c1 c3*s2*c1+s3*s1 0];
       [s3*c2   s3*s2*s1+c3*c1 s3*s2*c1-c3*s1 0];
       [-s2     c2*s1  c2*c1                  0]
       [0       0      0                      1]]; 
   
end

% ZYX rotztion
% [m] = [[c1*c2  -c2*s1  s2 0];
%        [c3*s1+s2*s3*c1  c3*c1-s3*s2*s1 -c2*s3 0];
%        [s3*s1-s2*c3*c1  c3*s1*s2+c1*s3 c2*c3  0];
%        [0               0              0      1]]; 

% YXZ rotztion
% [m] = [[c3*c1-s3*s2*s1  -c2*s3  c3*s1+s2*s3*c1 0];
%        [c3*s1*s2+c1*s3  c2*c3 s3*s1-s2*c3*c1 0];
%        [-c2*s1  s2 c2*c1  0];
%        [0       0  0      1]]; 

% ZXY rotztion
% m=[[c3*c1+s3*s2*s1  c1*s3*s2-c3*s1  c2*s3 0];
%    [c2*s1  c2*c1 -s2 0];
%    [c3*s1*s2-c1*s3  c1*c3*s2+s3*s1 c2*c3  0];
%    [0               0              0      1]]; 
