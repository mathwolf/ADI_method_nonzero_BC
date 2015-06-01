function IDM_CN_example_1()
% Summary of this function goes here
%   Detailed explanation goes here

% Use a uniform spatial grid of 0.01 in both x and y
N = 101;
%N = 5;
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


% Matrices used in ADI method
% This Bh d/n have the same constants as assignment 1
Bh = (1/h^2) * spdiags([-ones(N-2,1), 2*ones(N-2,1), -ones(N-2,1)], ...
                        [-1,0,1], N-2, N-2);
R = chol(speye(N-2) + (a*tau/2)*Bh);

% Matrix used to store alternating steps
U_twiddles = sparse(N-2,N-2);

% Vector for forcing term
load_vector = sparse(N-2,1);

% Step through the scheme
for m = 1:M
    % Create an array holding value of the forcing term, used on 
    % next step
    for i = 1:N-2
       for j = 1:N-2
          load_vector(i,j) = (1/2)*(f(h*i, h*j, tau*(m-1)) + ...
              f(h*i, h*j, tau*m)); 
       end
    end
    %{
    if m == 1
        disp(full(load_vector));
    end
    %}
    % Find the next U_twiddles 
    % Need to solve system once for each i = 1,...,N-1
    
    for i = 1:N-2  
        if i == 1
           b = (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
               (a*tau)/(2*h^2) * U(:,i+1) + ...
               (tau/2)*load_vector(:,i);
        elseif i == N-2
           b = (a*tau)/(2*h^2) * U(:,i-1) + ...
               (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
               (tau/2)*load_vector(:,i);
        else
           b = (a*tau)/(2*h^2) * U(:,i-1) + ...
               (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
               (a*tau)/(2*h^2) * U(:,i+1) + ...
               (tau/2)*load_vector(:,i);
        end
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
        U_twiddles(:,i) = R\b;
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
    % Now find next U
    for i = 1:N-2
       b = 2*U_twiddles(i,:)' - (speye(N-2) - (a*tau/2)*Bh) *  U(i,:)';
       %{
       if m == 1
          disp(full(b));
       end
       %}
       b = R'\b;
       U(i,:) = (R\b)';       
       %{
       if m == 1
          disp(full(U(i,:))); 
       end
       %}
    end
    %{
    if m ==1
       disp(full(U)); 
    end
    %}
    
    % Store data for plotting on selected steps
    
    %{
    if m == 1
        plot_data(2:N-1,2:N-1,2) = U(:,:);
    elseif m == 5
        plot_data(2:N-1,2:N-1,3) = U(:,:);
    %}
    
    
    if m == 1
        for i = 2:N-1
           for j = 2:N-1
              plot_data(i,j,2) = U(i-1,j-1) - u(h*(i-1), h*(j-1), tau*m); 
           end
        end
    elseif m == 5
        for i = 2:N-1
           for j = 2:N-1
              plot_data(i,j,3) = U(i-1,j-1) - u(h*(i-1), h*(j-1), tau*m); 
           end
        end 
    elseif m == 15
        for i = 2:N-1
           for j = 2:N-1
              plot_data(i,j,4) = U(i-1,j-1) - u(h*(i-1), h*(j-1), tau*m); 
           end
        end 
    elseif m ==40
        for i = 2:N-1
           for j = 2:N-1
              plot_data(i,j,5) = U(i-1,j-1) - u(h*(i-1), h*(j-1), tau*m); 
           end
        end 
        
        
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
%title('Approximation');

figure(2)
surf(X,Y,plot_data(:,:,2));
colormap winter;
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

