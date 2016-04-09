function x = psi2(y)

% Constants used to switch between different domains for testing.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;
global domain

% Upper boundary of test domain wrt y
if domain == SPHERE
    x = sqrt(1. - y^2);
elseif domain == ELLIPSOID
    x = sqrt(1. - (2.*y)^2);
elseif domain == ELLIPSE_CYLINDER
    x = sqrt(1. - (2.*y)^2);
elseif domain == CUBE
    x = 1;
elseif domain == OCTAHEDRON
    if (-1 <= y) && (y <= 0)
        x =  y + 1;
    elseif (0 <= y) && (y <= 1)
        x = - y + 1;
    else
        x = 0;
    end
elseif domain == DIAMOND_2
    if (-1 <= y) && (y <= 0)
        x =  0.5*y + 0.5;
    elseif (0 <= y) && (y <= 0.5)
        x = - y + 0.5;
    else
        x = 0;
    end
elseif domain == ELL
    if (-1 <= x) && (x <= 0)
        y = 0;
    elseif (0 < x) && (x <= 1)
        y = 1;
    else
        y = 0;
    end
elseif domain == RECTANGLE
    x = 1.;
elseif domain == CIRCLE_CYLINDER
    x = sqrt(1 - y^2);
end
end

