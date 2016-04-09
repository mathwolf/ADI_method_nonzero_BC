function y = theta2(x,z)

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

if domain == SPHERE
    y = sqrt(1. - x^2 - z^2);
elseif domain == CIRCLE_CYLINDER
    y = sqrt(1 - x^2);
elseif domain == ELLIPSE_CYLINDER
    y = 0.5 * sqrt(1. - x^2);
elseif domain == ELLIPSOID
    y = 0.5 * sqrt(1. - x^2 - 16.*z^2);
elseif domain == CUBE
    y = 1;
elseif domain == OCTAHEDRON
    if z >= 0
        if x >= 0
            y = - z - x + 1;
        else
            y = - z + x + 1;
        end
    else
        if x >= 0
            y = z - x + 1;
        else
            y = z + x + 1;
        end
    end
elseif domain == ELL
    if z <= 0
        if x <= 0
            y = 1;
        else
            y = 1;
        end
    else
        if x <= 0
            y = 0;
        else
            y = 1;
        end
    end
elseif domain == RECTANGLE
    y = 0.5;
end

end

