function y = phi2(x)

% Constants used to switch between different domains for testing.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
global domain

% Upper boundary of test domain wrt x
if domain == CIRCLE
    y = sqrt(1. - x^2);
elseif domain == ELLIPSE
    y = 0.5 * sqrt(1. - x^2);
elseif domain == DIAMOND
    if (-1 <= x) && (x <= 0)
        y = x + 1;
    elseif (0 <= x) && (x <= 1)
        y = - x + 1;
    else
        y = 0;
    end
elseif domain == DIAMOND_2
    if (-1 <= x) && (x <= 0)
        y = 0.5*x + 0.5;
    elseif (0 <= x) && (x <= 0.5)
        y = -x + 0.5;
    else
        y = 0;
    end
elseif domain == ELL
    if (-1 <= x) && (x <= 1)
        y = 1;
    else
        y = 0;
    end
elseif domain == RECTANGLE
    y = 0.5;
end
end

