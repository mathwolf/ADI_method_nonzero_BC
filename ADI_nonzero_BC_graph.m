function ADI_nonzero_BC_graph()
% Summary of this function goes here
%   Detailed explanation goes here

% Use a uniform spatial grid of 0.01 in both x and y
%N = 101;
N = 5;
h = 1./(N-1);

% Use a temporal spacing of the same size, go from 0 to 1
tau = 0.01;
M = 101;

% Parameter for material
a = 1.;

% Array of points used for plots
plot_data = zeros(N,N,5);

% Set up array of gridpoints and fill in with initial conditions
% Use initial conditions for first step
U = sparse(N,N);
for i = 1:N
   for j = 1:N
      U(i,j) = g1(h*(i-1), h*(j-1)); 
   end
end

%disp(full(U));             

% Store initial condition for plotting
for i = 1:N
   for j = 1:N
      plot_data(i,j,1) = U(i,j); 
%      plot_data(i,j,1) = U(i-1,j-1); 
   end
end


% Matrix used on left hand side of ADI method
% Use Cholesky factorization to solve system efficiently
Bh = (1/h^2) * spdiags([-ones(N-2,1), 2*ones(N-2,1), -ones(N-2,1)], ...
                        [-1,0,1], N-2, N-2);
R = chol(speye(N-2) + (a*tau/2)*Bh);

% Matrix used for alternating direction
U_twiddles = sparse(N,N);

% Array for forcing term
load_vector = sparse(N,N);

boundary_term = sparse(N-2,1);


% Evaluate each timestep
for m = 1:M
    % Create an array for the value of the forcing term at each 
    % gridpoint at one half of a timestep
    for i = 1:N
       for j = 1:N
          load_vector(i,j) = (1/2)*(f(h*(i-1), h*(j-1), tau*(m-1)) + ...
              f(h*(i-1), h*(j-1), tau*m)); 
       end
    end
    %{
    if m == 1
        disp(full(load_vector));
    end
    %}
    
    % Find the boundary values for the alternating direction.  Only the
    % values along the x boundary are needed by the formulas.
    for j = 2:N-1
       U_twiddles(1,j) = -a/(4*tau*h^2) * g2(0, h*(j-2), tau*m) + ...
           (1/2 + a/(2*tau*h^2)) * g2(0, h*(j-1), tau*m) + ...
           -a/(4*tau*h^2) * g2(0, h*j, tau*m) + ...
           -a/(4*tau*h^2) * g2(0, h*(j-2), tau*(m-1)) + ...
           (1/2 + a/(2*tau*h^2)) * g2(0, h*(j-1), tau*(m-1)) + ...
           -a/(4*tau*h^2) * g2(0, h*j, tau*(m-1));
       U_twiddles(N,j) = -a/(4*tau*h^2) * g2(1, h*(j-2), tau*m) + ...
           (1/2 + a/(2*tau*h^2)) * g2(1, h*(j-1), tau*m) + ...
           -a/(4*tau*h^2) * g2(1, h*j, tau*m) + ...
           -a/(4*tau*h^2) * g2(1, h*(j-2), tau*(m-1)) + ...
           (1/2 + a/(2*tau*h^2)) * g2(1, h*(j-1), tau*(m-1)) + ...
           -a/(4*tau*h^2) * g2(1, h*j, tau*(m-1));       
    end
    
    % Find the alternate matrix U_twiddles 
    % Need to solve system once for each i = 2,...,N-1
    for i = 2:N-1  
        
        boundary_term(1,1) = U_twiddles(1,i);
        boundary_term(N-2,1) = U_twiddles(N,i);
        b = (a*tau)/(2*h^2) * U(2:N-1,i-1) + ...
            (1 - 2*(a*tau)/(2*h^2)) * U(2:N-1,i) + ...
            (a*tau)/(2*h^2) * U(2:N-1,i+1) + ...
            (tau/2)*load_vector(2:N-1,i) + ...
            (a*tau)/(2*h^2) * boundary_term(:,1);

        %{
        b = (speye(N-2) - (a*tau/2)*Bh) * U(i,:)' + ...
            (tau/2)*load_vector(i,:)';
        %}
        %{
        if m == 1
           disp(full((speye(N-2) - (a*tau/2)*Bh) * U(i,:)'));
           disp(full(b)); 
        end
            %}
        b = R'\b;
        %{
        if m == 1
           disp(full(b)); 
        end
        %}
        U_twiddles(2:N-1,i) = R\b;
        %{
        if m ==1
           disp(full(U_twiddles(:,i))); 
        end
        %}
    end
    %{
    if m == 1
       disp(full(U_twiddles)); 
    end
    %}
    
    % Now find the interior points for next U
    % Solve system one for each j = 2, ... , N-1
    for j = 2:N-1
        
       boundary_term(1,1) = g2(h*(j-1), 0, tau*m) - U(j,1);
       boundary_term(N-2,1) = g2(h*(j-1), 1, tau*m) - U(j,N);
       b = 2*U_twiddles(i,2:N-1)' - ...
           (speye(N-2) - (a*tau/2)*Bh) *  U(i,2:N-1)' + ...
           (a*tau)/(2*h^2) * boundary_term(:,1);
       %{
       if m == 1
          disp(full(b));
       end
       %}
       b = R'\b;
       U(i,2:N-1) = (R\b)';       
       %{
       if m == 1
          disp(full(U(i,:))); 
       end
       %}
    end
    
    % Update boundary points
    for i = 1:N
       U(1,i) = g2(0, h*(i-1), tau*m); 
       U(N,i) = g2(1, h*(i-1), tau*m);
       U(i,1) = g2(h*(i-1), 0, tau*m);
       U(i,N) = g2(h*(i-1), 1, tau*m);
    end
    
    %{
    if m ==1
       disp(full(U)); 
    end
    %}
    
    % Store data for plotting on selected steps
      
    
    if m == 1
        for i = 1:N
           for j = 1:N
              plot_data(i,j,2) = U(i,j); 
           end
        end
    elseif m == 5
        for i = 1:N
           for j = 1:N
              plot_data(i,j,3) = U(i,j); 
           end
        end 
    elseif m == 30
        for i = 1:N
           for j = 1:N
              plot_data(i,j,4) = U(i,j); 
           end
        end 
    elseif m == 70
        for i = 1:N
           for j = 1:N
               %{
              plot_data(i,j,5) = U(i,j) - ...
                  u(h*(i-1), h*(j-1), tau*(m-1));
                %}
               plot_data(i,j,5) = U(i,j);
           end
        end 
        
  %{  
    if m == 1
        for i = 1:N
           for j = 1:N
              plot_data(i,j,2) = U(i,j) - ...
                  u(h*(i-1), h*(j-1), tau*(m-1)); 
           end
        end
    elseif m == 5
        for i = 1:N
           for j = 1:N
              plot_data(i,j,3) = U(i,j) - ...
                  u(h*(i-1), h*(j-1), tau*(m-1)); 
           end
        end 
    elseif m == 30
        for i = 1:N
           for j = 1:N
              plot_data(i,j,4) = U(i,j) - ...
                  u(h*(i-1), h*(j-1), tau*(m-1)); 
           end
        end 
    elseif m == 70
        for i = 1:N
           for j = 1:N
               %{
              plot_data(i,j,5) = U(i,j) - ...
                  u(h*(i-1), h*(j-1), tau*(m-1));
                %}
               plot_data(i,j,5) = U(i,j);
           end
        end 
%}
        
        
        %{
    elseif m == 15
        plot_data(2:N-1,2:N-1,4) = U(:,:);
    elseif m == 40
        plot_data(2:N-1,2:N-1,5) = U(:,:);
        %}
    end
end

% Display exact solution at first timestep
%{
V = sparse(N-2,N-2);
for i = 1:N-2
    for j = 1:N-2
        V(i,j) = u(h*i, h*j, tau/2);
    end
end
U = plot_data(2:N-1,2:N-1,2);
disp(full(U));
disp(full(V));
disp(full(U-V));
%}

% Plot of exact solution at first timestep
%{
for i = 2:N-1
   for j = 2:N-1
      plot_data(i,j,2) = u(h*(i-1), h*(j-1),0.01); 
   end
end
%}

% Create plot
X = linspace(0., 1., N);
Y = linspace(0., 1., N);

% Make plots of error at each timestep for debugging

figure(1)
surf(X,Y,plot_data(:,:,1));
colormap winter;
xlabel('x');
ylabel('y');
%title('Approximation');

figure(2)
surf(X,Y,plot_data(:,:,2));
colormap winter;
xlabel('x');
ylabel('y');

%title('Exact solution');


figure(3)
surf(X,Y,plot_data(:,:,3));
colormap winter;

figure(4)
surf(X,Y,plot_data(:,:,4));
colormap winter;

figure(5)
surf(X,Y,plot_data(:,:,5));
colormap winter;

end

