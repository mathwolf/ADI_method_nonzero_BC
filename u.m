function z = u(x,y,t)

% Constants used to select different testing parameters
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
EXPONENT_3 = 3;
TRIG = 4;

% Exact solution used as test function
global test_solution

if test_solution == EXPONENT_0
    z = x*(1-x)*y*(1-y)*exp(x+y+t);
elseif test_solution == EXPONENT_1
    z = exp(x+y+t);
elseif test_solution == EXPONENT_2
    z = exp(x + 2*y + 3*t);
elseif test_solution == EXPONENT_3
    z = exp(x*y*t);
elseif test_solution == TRIG
    z = 2 * sin(10*(x^2+y^2+t^2));
else
    z = 0;
end
end

