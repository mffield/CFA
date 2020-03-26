function Q = M2Q(R)
% This function converts a rotation matrix into quaternion (euler
% parameter) form.
% Input: Rotation matrix 3x3
% Output: quaternion 1x4

% Written by Matthew Field

if length(R)==4
    R = R(1:3,1:3);
end

Q = zeros(1,4);
Q(4) = (1/2)*(1 + R(1,1) + R(2,2) + R(3,3))^(0.5);
Q(1) = (1/2)*(1 + R(1,1) - R(2,2) - R(3,3))^(0.5);
Q(2) = (1/2)*(1 - R(1,1) + R(2,2) - R(3,3))^(0.5);
Q(3) = (1/2)*(1 - R(1,1) - R(2,2) + R(3,3))^(0.5);
[e,I]=max(Q);
switch I
    case 1
        Q(2) = (R(1,2) + R(2,1))/(4*Q(1));
        Q(3) = (R(1,3) + R(3,1))/(4*Q(1));
        Q(4) = (R(3,2) - R(2,3))/(4*Q(1));
    case 2
        Q(1) = (R(1,2) + R(2,1))/(4*Q(2));
        Q(3) = (R(3,2) + R(2,3))/(4*Q(2));
        Q(4) = (R(1,3) - R(3,1))/(4*Q(2));
    case 3
        Q(1) = (R(1,3) + R(3,1))/(4*Q(3));
        Q(2) = (R(3,2) + R(2,3))/(4*Q(3));
        Q(4) = (R(2,1) - R(1,2))/(4*Q(3));
    case 4
        Q(1) = (R(3,2) - R(2,3))/(4*Q(4));
        Q(2) = (R(1,3) - R(3,1))/(4*Q(4));
        Q(3) = (R(2,1) - R(1,2))/(4*Q(4));
end
Q = [Q(4) Q(1:3)];
end

