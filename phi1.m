function y = phi1(x)

% Constants used to switch between different domains for testing.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;

CIRCLE_CYLINDER = 8;
ELLIPSE_CYLINDER = 9;
global domain

% Lower boundary of test domain wrt x
if domain == SPHERE
    y = - sqrt(1. - x^2);
elseif domain == ELLIPSOID
    y = - 0.5 * sqrt(1. - x^2);
elseif domain == ELLIPSE_CYLINDER
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
    if (-1 <= x) && (x <= 0)
        y = -1;
    elseif (0 < x) && (x <= 1)
        y = 0;
    else
        y = 0;
    end
elseif domain == RECTANGLE
    y = -0.5;
elseif domain == CIRCLE_CYLINDER
    y = -sqrt(1 - x^2);
end
end

