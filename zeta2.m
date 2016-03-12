function z = zeta2(x,y)

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
    z = sqrt(1. - x^2 - y^2);
elseif domain == ELLIPSOID
    z = 0.25 * sqrt(1. - x^2 - 4.*y^2);
elseif domain == CUBE
    z = 1;
elseif domain == OCTAHEDRON
    if x >= 0
        if y >= 0
            z = - x - y + 1;
        else
            z = - x + y + 1;
        end
    else
        if y >= 0
            z = x - y + 1;
        else
            z = x + y + 1;
        end
    end
elseif domain == ELL
    if x <= 0
        if y <= 0
           z = 0; 
        else
           z =0;
        end
    else
        if y <= 0
           z = 1; 
        else
           z = 1;
        end
    end
    
elseif domain == RECTANGLE
    z = 0.25;
end

end

