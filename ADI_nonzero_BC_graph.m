function ADI_nonzero_BC_graph()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions. This program tests the method using the function
% e^(x+y+t) using spatial and temporal grids with spacing 0.01.  The
% exact solution, approximate solution, and error are displayed
% graphically.

% Use a uniform spatial grid of 0.01 in both x and y
N = 100;
h = 1./N;

% Use a temporal spacing of the same size, go from 0 to 1
tau = 0.01;
M = 101;

% Array of points used for plots
plot_data = zeros(N+1,N+1,6);

% Set up gridpoints and fill in with initial conditions
U = zeros(N+1,N+1);
for i = 0:N
   for j = 0:N
      U(i+1, j+1) = g1(h*i, h*j); 
   end
end

% Create tridiagonal matrix used on left hand side of ADI method
tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
A = speye(N-1) + (tau/(2*h^2)) * tridiagonal;

% Matrix V used on RHS of matrix formulas, interior points only
V = zeros(N-1,N-1);

% Array of values for forcing term
load_array = zeros(N+1,N+1);

% Evaluate each timestep
for m = 1:M
    % Create an array for the forcing term at each 
    % gridpoint at one half of a timestep
    for i = 1:N-1
       for j = 1:N-1
          load_array(i,j) = (1/2)*(f(h*i, h*j, tau*(m-1)) + ...
              f(h*i, h*j, tau*m)); 
       end
    end
    
    % Create an array of values for V at each interior gridpoint
    for i = 1:N-1
       for j = 1:N-1
          V(i,j) = U(i+1,j+1) + ...
              (tau / (2*h^2)) * (U(i+1,j) - 2*U(i+1,j+1) + U(i+1,j+2));
       end
    end
    
   
    % Find the approximation at the half timestep
    % Need to solve system once for each j = 2,...,N
    for j = 1:N-1  
        % Find boundary values for first and last point, no perturbation
        % term
        btf = (1/2) * (g2(0, h*j, tau*m) + g2(0, h*j, tau*(m-1)));
        btl = (1/2) * (g2(1, h*j, tau*m) + g2(1, h*j, tau*(m-1)));
        % Find rhs of vector equation
        b = V(1:N-1,j) + (tau/2) * load_array(1:N-1,j);
        % Add effect of the boundary points
        b(1) = b(1) + (tau/(2*h^2)) * btf;
        b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
        % Solve for approximate values at half timestep
        U(2:N,j+1) = A\b;
    end
    
    % Now find the interior points for next whole timestep
    % Solve system one for each i = 2, ... , N
    for i = 1:N-1        
       % Find boundary values for first and last point
       btf = g2(h*i, 0, tau*m);
       btl = g2(h*i, h*N, tau*m);
       % Get the rhs of the vector equation
       b = 2*U(i+1,2:N)' - V(i,1:N-1)';
       % Add effect of the boundary points
       b(1) = b(1) + (tau/(2*h^2)) * btf;
       b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
       % Solve vector equation for the next full timestep
       U(i+1,2:N) = (A\b)';       
    end
    
    % Update boundary points
    for i = 0:N
       U(1,i+1) = g2(0, h*i, tau*m); 
       U(N+1,i+1) = g2(h*N, h*i, tau*m);
       U(i+1,1) = g2(h*i, 0, tau*m);
       U(i+1,N+1) = g2(h*i, h*N, tau*m);
    end
    
    % Store data for plotting on selected steps    
    if m == 1
        for i = 1:N+1
           for j = 1:N+1
              plot_data(i,j,1) = U(i,j);
              plot_data(i,j,2) = u(h*(i-1), h*(j-1), tau*m);
              plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
           end
        end
    elseif m == M
        for i = 1:N+1
           for j = 1:N+1
              plot_data(i,j,4) = U(i,j);
              plot_data(i,j,5) = u(h*(i-1), h*(j-1), tau*m);
              plot_data(i,j,6) = plot_data(i,j,1) - plot_data(i,j,2);
           end
        end 
    end
end

% Create plots
X = linspace(0., 1., N+1);
Y = linspace(0., 1., N+1);

figure
subplot(2,3,1)
surf(X,Y,plot_data(:,:,1))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after one timestep')

subplot(2,3,2)
surf(X,Y,plot_data(:,:,2))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after one timestep')

subplot(2,3,3)
surf(X,Y,plot_data(:,:,3));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after one timestep')

subplot(2,3,4)
surf(X,Y,plot_data(:,:,4))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after final timestep')

subplot(2,3,5)
surf(X,Y,plot_data(:,:,5))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after final timestep')

subplot(2,3,6)
surf(X,Y,plot_data(:,:,6));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after final timestep')


end

