function z = f(x,y,t)
% Forcing term for test function

EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;

global solution

if solution == TRIG
    z = 10 * (1024*x^2 + 64*y^2) * cos(16*x^2 + 4*y^2 + t) + ...
        390 * sin(16*x^2 + 4*y^2 + t);
elseif solution == EXPONENT_1
    z = - exp(x+y+t);
elseif solution == EXPONENT_2
    z = - 12 * exp(3*x + 2*y + t);
elseif solution == POLY
    z = 10 - 240*x + 96*y;
end

end

