function z = f(x,y,t)
% Forcing term for test function

z = -3 * cos(x^2 + y^2 + t) + (4*x^2 + 4*y^2) * sin(x^2 + y^2 + t);

end

