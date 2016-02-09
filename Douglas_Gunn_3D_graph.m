function Douglas_Gunn_3D_graph()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions in 3D.

FALSE = 0;
TRUE = 1;

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution
test_solution = EXPONENT_1; % exp test function with 0 BC

% Turn on or off perturbation on half step boundary conditions
perturbation = TRUE;

% Use a uniform spatial grid of 0.01 in x, y and z
N = 20;
h = 1./N;

% Use a temporal spacing of the same size, go from 0 to 1
tau = h;
M = N+1;

% Array of points used for plots
plot_data = zeros(N+1,N+1,6);

% Set up gridpoints and fill in with initial conditions
U = zeros(N+1,N+1,N+1);
for i = 0:N
   for j = 0:N
       for k = 0:N
            U(i+1, j+1, k+1) = g3(h*i, h*j, h*k); 
       end
   end
end

% Create tridiagonal matrix used on left hand side of ADI method
tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
A = speye(N-1) + (tau/(2*h^2)) * tridiagonal;

% Create tridiagonal matrix used on right hand side of method
D = speye(N-1) - (tau/(2*h^2)) * tridiagonal;

% Array of gridpoints V used on RHS of matrix equations
V = zeros(N+1,N+1,N+1);

% Array of values for forcing term, interior points only
load_array = zeros(N-1,N-1,N-1);

% Evaluate each timestep
for m = 1:M
    % Create an array for the forcing term at each interior
    % gridpoint at one half of a timestep
    for i = 1:N-1
       for j = 1:N-1
           for k = 1:N-1
            load_array(i,j,k) = (1/2)*(f3(h*i, h*j, h*k, tau*(m-1)) + ...
              f3(h*i, h*j, h*k, tau*m)); 
           end
       end
    end
    
    % Find only the boundary values of V (the difference array) at
    % the next full timestep.  We will use these values to get the
    % BCs on the subsequent steps of the method.
    for i = 0:N
       for j = 0:N
          V(1,i+1,j+1) = g4(0, h*i, h*j, tau*m) ...
              - g4(0, h*i, h*j, tau*(m-1));
          V(N+1, i+1, j+1) =  g4(h*N, h*i, h*j, tau*m)...
              - g4(h*N, h*i, h*j, tau*(m-1));
          V(i+1, 1, j+1) = g4(h*i, 0, h*j, tau*m) ...
              - g4(h*i, 0, h*j, tau*(m-1));
          V(i+1, N+1, j+1) = g4(h*i, h*N, h*j, tau*m) ...
              - g4(h*i, h*N, h*j, tau*(m-1));
          V(i+1,j+1,1) = g4(h*i, h*j, 0, tau*m) ...
              - g4(h*i, h*j, 0, tau*(m-1));
          V(i+1, j+1, N+1) = g4(h*i, h*j, h*N, tau*m) ...
              - g4(h*i, h*j, h*N, tau*(m-1));
       end
    end    
       
    % Find the approximation for V at the one-third timestep
    % Need to solve system once for each j = 2,...,N
    % and k = 2,...,N
    for j = 2:N
        for k = 2:N
            % Find boundary values for first and last point, no perturbation
            % term
            if perturbation == TRUE            
                btf = V(1,j,k) - (tau/(2*h^2)) * ...
                    (V(1,j,k-1) - 2 * V(1,j,k) + V(1,j,k+1)) ...
                    - (tau/(2*h^2)) * ...
                    (V(1,j-1,k) - 2 * V(1,j,k) + V(1,j+1,k)) ...
                    +(tau^2/(4*h^4)) * ...
                    (V(1,j-1,k-1) - 2 * V(1,j,k-1) + V(1,j+1,k-1) ...
                    - 2 * V(1,j-1,k) + 4 * V(1,j,k) - 2 * V(1,j+1,k) ...
                    + V(1,j-1,k+1) - 2 * V(1,j,k+1) + V(1,j+1,k+1));
                btl = V(N+1,j,k) - (tau/(2*h^2)) * ...
                    (V(N+1,j,k-1) - 2 * V(N+1,j,k) + V(N+1,j,k+1)) ...
                    - (tau/(2*h^2)) * ...
                    (V(N+1,j-1,k) - 2 * V(N+1,j,k) + V(N+1,j+1,k)) ...
                    +(tau^2/(4*h^4)) * ...
                    (V(N+1,j-1,k-1) - 2 * V(N+1,j,k-1) + V(N+1,j+1,k-1) ...
                    - 2 * V(N+1,j-1,k) + 4 * V(N+1,j,k) - 2 * V(N+1,j+1,k) ...
                    + V(N+1,j-1,k+1) - 2 * V(N+1,j,k+1) + V(N+1,j+1,k+1));
            else
                btf = V(1,j,k);
                btl = V(N+1,j,k);

            end
            
            % Find rhs of vector equation
            % First, create a vector describing the three terms that 
            % contain a difference quotients.
            b = zeros(N-1,1);
            for i = 2:N
                % Add effects of x derivative
                b(i-1) = (tau / (h^2)) * ...
                    (U(i-1,j,k) - 2*U(i,j,k) + U(i+1,j,k));
                % Add effects of y derivative
                b(i-1) = b(i-1) + (tau / (h^2)) * ...
                    (U(i,j-1,k) - 2*U(i,j,k) + U(i,j+1,k));
                % Add effects of z derviative
                b(i-1) = b(i-1) + (tau / (h^2)) * ...
                    (U(i,j,k-1) - 2*U(i,j,k) + U(i,j,k+1));
            end
            % Add the effect of the forcing term
            b = b + tau * load_array(1:N-1,j-1,k-1);
            % Add effect of the boundary points
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve for approximate values of V at one-third timestep
            V(2:N,j,k) = A\b;
        end
    end
    
    % Store plot data for debugging
%     for i = 1:N+1
%         for j = 1:N+1
%             plot_data(i,j,1) = U(i,j,3);
%             plot_data(i,j,2) = u3(h*(i-1),h*(j-1),0.5,tau*(m-(1./2.)));
%             plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
%         end
%     end
    
    % Find the approximation for V at the two-thirds timestep
    % Need to solve system once for each i = 2,...,N
    % and k = 2,...,N
    for i = 2:N
        for k = 2:N
            % Find boundary values for first and last point, no perturbation
            % term
            if perturbation == TRUE            
                btf = V(i,1,k) - (tau/(2*h^2)) * ...
                    (V(i,1,k-1) - 2 * V(i,1,k) + V(i,1,k+1));
                btl = V(i,N+1,k) - (tau/(2*h^2)) * ...
                    (V(i,N+1,k-1) - 2 * V(i,N+1,k) + V(i,N+1,k+1));
            else
                btf = V(i,1,k);
                btl = V(i,N+1,k);
            end
            
            % Find rhs of vector equation
            b = V(i,2:N,k)';
            % Add effect of the boundary points
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve for approximate values of V at two-thirds timestep
            V(i,2:N,k) = (A\b)';
        end
    end
    
    % Store plot data for debugging
%     for i = 1:N+1
%         for j = 1:N+1
%             plot_data(i,j,4) = U(i,j,3);
%             plot_data(i,j,5) = u3(h*(i-1),h*(j-1),0.5,tau*(m-(1./3.)));
%             plot_data(i,j,6) = plot_data(i,j,4) - plot_data(i,j,5);
%         end
%     end
        
    % Now find the approximation for V at next whole timestep
    % Solve system one for each fixed i and j
    for i = 2:N
        for j = 2:N
            btf = V(i,j,1);
            btl = V(i,j,N+1);
            
            % Get the rhs of the vector equation
            b = squeeze(V(i,j,2:N));
            % Add effect of the boundary points
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve vector equation for the next full timestep
            V(i,j,2:N) = permute(A\b, [3 2 1]);       
        end
    end
    
    % Get the next whole timestep for U, the approximation
    % to the actual solution of the PDE
    U = U + V;
    
    % Update boundary points
    for i = 0:N
       for j = 0:N
          U(1,i+1,j+1) = g4(0, h*i, h*j, tau*m);
          U(N+1, i+1, j+1) = g4(h*N, h*i, h*j, tau*m);
          U(i+1, 1, j+1) = g4(h*i, 0, h*j, tau*m);
          U(i+1, N+1, j+1) = g4(h*i, h*N, h*j, tau*m);
          U(i+1,j+1,1) = g4(h*i, h*j, 0, tau*m);
          U(i+1, j+1, N+1) = g4(h*i, h*j, h*N, tau*m);
       end
    end
    % Store plot data for debugging
%     for i = 1:N+1
%         for j = 1:N+1
%             plot_data(i,j,7) = U(i,j,3);
%             plot_data(i,j,8) = u3(h*(i-1),h*(j-1),0.5,tau*m);
%             plot_data(i,j,9) = plot_data(i,j,7) - plot_data(i,j,8);
%         end
%     end
        
    
    % Store data for plotting on selected steps    
    if m == 1
        for i = 1:N+1
           for j = 1:N+1
              plot_data(i,j,1) = U(i,j,2);
              plot_data(i,j,2) = u3(h*(i-1), h*(j-1), h, tau*m);
              plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
           end
        end
    elseif m == M
        for i = 1:N+1
           for j = 1:N+1
              plot_data(i,j,4) = U(i,j,2);
              plot_data(i,j,5) = u3(h*(i-1), h*(j-1), h, tau*m);
              plot_data(i,j,6) = plot_data(i,j,4) - plot_data(i,j,5);
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

% subplot(3,3,7)
% surf(X,Y,plot_data(:,:,7))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Approximate solution after final timestep')
% 
% subplot(3,3,8)
% surf(X,Y,plot_data(:,:,8))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Exact solution after final timestep')
% 
% subplot(3,3,9)
% surf(X,Y,plot_data(:,:,9));
% colormap winter;
% xlabel('x')
% ylabel('y')
% title('Error after final timestep')


end

