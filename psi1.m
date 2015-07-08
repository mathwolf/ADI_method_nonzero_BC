function x = psi1(y, domain)

% Constants used to switch between different domains for testing.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;

% Lower boundary of test domain wrt y
if domain == CIRCLE
    x = - sqrt(1. - y^2);
elseif domain == ELLIPSE
    x = - sqrt(1. - (2.*y)^2);
elseif domain == DIAMOND
    if (-1 <= y) && (y <= 0)
        x =  - y - 1;
    elseif (0 <= y) && (y <= 1)
        x = y - 1;
    else
        x = 0;
    end
elseif domain == ELL
    if (-1 <= y) && (y <= 0)
        x = 0;
    elseif (0 <= y) && (y <= 1)
        x = -1;
    else
        x = 0;
end
end

