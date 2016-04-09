function x = eta2(y,z)

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
    x = sqrt(1. - y^2 - z^2);
elseif domain == CIRCLE_CYLINDER
    x = sqrt(1 - y^2);
elseif domain == ELLIPSE_CYLINDER
    x = sqrt(1. - (2.*y)^2);
elseif domain == ELLIPSOID
    x = sqrt(1. - 4.*y^2 - 16.*z^2);
elseif domain == CUBE
    x = 1.;
elseif domain == OCTAHEDRON
    if y >= 0
        if z >= 0
            x = - y - z + 1;
        else
            x = - y + z + 1;
        end
    else
        if z >= 0
            x = y - z + 1;
        else
            x = y + z + 1;
        end
    end
elseif domain == ELL
    if y <= 0
        if z <= 0
            x = 0;
        else
            x = 0;
        end
    else
        if z <= 0
            x = 1;
        else
            x = 1;
        end
    end
elseif domain == RECTANGLE
    x = 1.;
end

end
