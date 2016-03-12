function x = psi1(y)

% Constants used to switch between different domains for testing.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;
global domain

% Lower boundary of test domain wrt y
if domain == SPHERE
    x = - sqrt(1. - y^2);
elseif domain == ELLIPSOID
    x = - sqrt(1. - (2.*y)^2);
elseif domain == CUBE
    x = 0;
elseif domain == OCTAHEDRON
    if (-1 <= y) && (y <= 0)
        x =  - y - 1;
    elseif (0 <= y) && (y <= 1)
        x = y - 1;
    else
        x = 0;
    end
elseif domain == DIAMOND_2
    if (-1 <= y) && (y <= 0)
        x =  - y - 1;
    elseif (0 <= y) && (y <= 0.5)
        x = 2*y - 1;
    else
        x = 0;
    end
elseif domain == ELL
    x = -1.;
elseif domain == RECTANGLE
    x = -1.;
end
end
