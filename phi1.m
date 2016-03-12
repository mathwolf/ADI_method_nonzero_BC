function y = phi1(x)

% Constants used to switch between different domains for testing.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;
global domain

% Lower boundary of test domain wrt x
if domain == SPHERE
    y = - sqrt(1. - x^2);
elseif domain == ELLIPSOID
    y = - 0.5 * sqrt(1. - x^2);
elseif domain == CUBE
    y = 0;
elseif domain == OCTAHEDRON
    if (-1 <= x) && (x <= 0)
        y = -x - 1;
    elseif (0 <= x) && (x <= 1)
        y = x - 1;
    else
        y = 0;
    end
elseif domain == DIAMOND_2
    if (-1 <= x) && (x <= 0)
        y = -x - 1;
    elseif (0 <= x) && (x <= 0.5)
        y = 2*x - 1;
    else
        y = 0;
    end
elseif domain == ELL
    y = -1;
elseif domain == RECTANGLE
    y = -0.5;
end
end

