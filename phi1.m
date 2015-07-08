function y = phi1(x, domain)

% Constants used to switch between different domains for testing.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;

% Lower boundary of test domain wrt x
if domain == CIRCLE
    y = - sqrt(1. - x^2);
elseif domain == ELLIPSE
    y = - sqrt((0.5)*(1. - x^2));
elseif domain == DIAMOND
    if (-1 <= x) && (x <= 0)
        y = -x - 1;
    elseif (0 <= x) && (x <= 1)
        y = x - 1;
    else
        y = 0;
    end
elseif domain == ELL
    if (-1 <= x) && (x <= 0)
        y = 0;
    elseif (0 <= x) && (x <= 1)
        y = - 1;
    else
        y = 0;
    end
end
end

