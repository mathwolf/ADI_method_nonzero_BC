function z = f(x,y,t)
% Constants used to select different testing parameters
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
EXPONENT_3 = 3;
TRIG = 4;

% Forcing term for test function
global test_solution

if test_solution == EXPONENT_0
    z = -x*y*(x*y+3*x+3*y-7)*exp(x+y+t);
elseif test_solution == EXPONENT_1
        z = -exp(x+y+t);
elseif test_solution == EXPONENT_2
    z = -2 * exp(x + 2*y + 3*t);
elseif test_solution == EXPONENT_3
    z = (x*y - x^2*t^2 - y^2*t^2) * exp(x*y*t);
elseif test_solution == TRIG
    z = 40*(t-2)*cos(10*(x^2+y^2+t^2)) + ...
        800*(x^2+y^2)*sin(10*(x^2+y^2+t^2));
else
    z = 0;
end
end

