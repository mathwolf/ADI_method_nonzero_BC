function y = theta1(x,z)

% Constants used to switch between different domains for testing.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;
global domain

if domain == SPHERE
    y = - sqrt(1. - x^2 - z^2);
elseif domain == ELLIPSOID
    y = - 0.5 * sqrt(1. - x^2 - 16.*z^2);
elseif domain == CUBE
    y = 0;
elseif domain == OCTAHEDRON
    if z >= 0
        if x >= 0
            y = z + x - 1;
        else
            y = z - x - 1;
        end
    else
        if x >= 0
            y = - z + x - 1;
        else
            y = - z - x - 1;
        end
    end
elseif domain == RECTANGLE
    y = -0.5;
end

end

