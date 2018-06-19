function [jA_x,jA_y,jA_z] = jA3(x)

xCenter = pi/2;
c = pi/10;
gauss = 1*exp(- (x-xCenter)^2 / (2*c^2) );

jA_x = gauss * complex(1,1);
jA_y = gauss * complex(1,1);
jA_z = gauss * complex(1,1);

end
