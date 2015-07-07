function z = u(x,y,t)
% Exact solution u(x,y,t) used as test function

EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;

global solution

if solution == TRIG
    z = 10 * cos(16*x^2 + 4*y^2 + t);
elseif solution == EXPONENT_1
    z = exp(x + y + t);
elseif solution == EXPONENT_2
    z = exp(3*x + 2*y + t);
elseif solution == POLY
    z = 40*x^3 - 60*x^2*y + 20*y^3 + 10*t;
end

end

