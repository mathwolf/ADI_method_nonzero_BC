function x = eta2(y,z)

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
    x = sqrt(1. - y^2 - z^2);
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
elseif domain == RECTANGLE
    x = 1.;
end

end
