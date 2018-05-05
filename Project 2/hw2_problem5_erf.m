% This script performs all testing needed to analyze the three ode solvers.

% The differential equation given in Problem 5
f = @(x,y) (2/sqrt(pi)) * exp(-x^2);

tols = [1e-1, 1e-4, 1e-7, 1e-10, 1e-13];
solvers = [23, 45, 113]; % to help format output in the loop

results = cell(3,5); % stores the results from each solver at each rel. tol.
errors  = cell(3,5); % stores the errors from each solver at each rel. tol.
steps   = zeros(3,5); % stores the number of steps from each solver at each rel. tol.
times   = zeros(3,5); % stores the times from each solver at each rel. tol.

% for each of the 5 rel. tols...
for j=1:5
    % Stores the results from each solver at each rel. tol.
    results{1,j} = ode23(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    results{2,j} = ode45(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    results{3,j} = ode113(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    
    % Stores the ode function calls as anonymous handles for timing
    temp_23  = @() ode23(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    temp_45  = @() ode45(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    temp_113 = @() ode113(f, [0, 2], erf(0), odeset('RelTol', tols(j)));
    
    % Store the times of each ode function call
    times(1,j) = timeit(temp_23);
    times(2,j) = timeit(temp_45);
    times(3,j) = timeit(temp_113);
    
    % Display the current rel. tol. to organize the ouput
    disp(['RelTol: ', num2str(tols(j))]);
    
    % for each of the 3 ode solvers...
    for i=1:3
        % Store the errors at each point for each solver
        errors{i,j} = abs(results{i,j}.y - erf(results{i,j}.x))';
        
        % Output the maximum error and endpoint error for each solver
        disp(['ode' num2str(solvers(i)) ' max error: '  num2str(max(errors{i,j}))]);
        disp(['ode' num2str(solvers(i)) ' end error: '  num2str(errors{i,j}(end))]);
        
        % Store the number of steps taken by each solver
        steps(i,j) = length(errors{i,j});
    end
    
    disp(' '); % To separate the outputs from each rel. tol. with a newline
    
    % Plots the errors at each step for each solver
    figure
    plot(results{1,j}.x, errors{1,j}, results{2,j}.x, errors{2,j}, results{3,j}.x, errors{3,j});
    title(['Errors of Solvers with a ', num2str(tols(j)), ' Relative Tolerance'])
    xlabel('0 < x < 2') % x-axis label
    ylabel('Local Errors of ODE Solvers') % y-axis label
    legend('ode23','ode45','ode113', 'Location','northwest')
end

% Produces a barchart of the steps taken for each solver at each rel. tol.
figure
bar(categorical(tols), steps.')
title('Steps Taken by Solvers with Various Relative Tolerances')
xlabel('Relative Tolerances')
ylabel('Number of Steps Taken') % y-axis label
legend('ode23','ode45','ode113', 'Location','northeast')

% Produces a barchart of the time taken for each solver at each rel. tol.
figure
bar(categorical(tols), times.')
title('Time Taken by Solvers with Various Relative Tolerances')
xlabel('Relative Tolerances')
ylabel('Time (in seconds)') % y-axis label
legend('ode23','ode45','ode113', 'Location','northeast')
