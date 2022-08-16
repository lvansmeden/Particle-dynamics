function Rrot = RotationMatrix(d_theta)
Rx = [[1 0 0];[0 cos(d_theta(1)) -sin(d_theta(1))];[0 sin(d_theta(1)) cos(d_theta(1))]];
Ry = [[cos(d_theta(2)) 0 sin(d_theta(2))];[0 1 0];[-sin(d_theta(2)) 0 cos(d_theta(2))]];
Rz = [[cos(d_theta(3)) -sin(d_theta(3)) 0];[sin(d_theta(3)) cos(d_theta(3)) 0];[0 0 1]];

Rrot = Rx*Ry*Rz;
end 


