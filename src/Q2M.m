function [RT, X] = Q2M(pos, q)   
%calculate generalised rotation matrices from position and quaternion
%vectors

% Input: position vector of length 3, quaternion vector of length 4 where q = [w,x,y,z]
% Output: generalised rotation matrix RT 4x4, and rotation matrix 3x3

% Written by Matthew Field
    
    w=q(1); x=q(2); y=q(3); z=q(4);
%     w=q(4); x=q(1); y=q(2); z=q(3);

    xx = x^2;    yy = y^2;    zz = z^2;
    wx = w*x;    xy = x*y;    wz = w*z;    wy = w*y;    yz = y*z; 
    xz = x*z;    
    RT=[1-2*(yy+zz), 2*(xy-wz), 2*(xz+wy), pos(1);
        2*(xy+wz), 1-2*(xx+zz), 2*(yz-wx), pos(2);
        2*(xz-wy), 2*(yz+wx), 1-2*(xx+yy), pos(3);
        0,0,0,1];
    
    
    X=RT(1:3,1:3);
end