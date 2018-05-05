% For Section 5.7 Problem 2(b)
f1 = @(t,y) sin(t) + exp(-t);
disp('Results for Section 5.7 Problem 2(b)');
adams_variable_pc(f1, 0, 1, 0, 1e-4, 0.2, 0.01);

disp(newline);

% For Section 5.7 Problem 7
f3 = @(t,y) (2.9 * 1e-2) * y - ((1.4 * 1e-7) * y^2);
disp('Results for Section 5.7 Problem 7');
adams_variable_pc(f3, 0, 5, 50976, 1e-1, 1, 0.05);


% For testing (Section 5.7 Example 1)
% f2 = @(t,y) y - t^2 + 1;
% adams_variable_pc(f2, 0, 2, 0.5, 1e-5, 0.2, 0.01);

% Adams Variable Step-Size Predictor-Corrector
% This function performs the Adams Variable Step-Size Predictor-Corrector
% algorithm similar to what is described in section 5.7.
function [i, t, w, h] = adams_variable_pc(f, a, b, y_a, TOL, hmax, hmin)
    %disp("Step 2");
    t = zeros(4,1); w = zeros(4,1);
    t(1) = a; w(1) = y_a; h = hmax; FLAG = 1; LAST = 0;
    disp([0 t(1) w(1)]);
    
    %disp("Step 3");
    [t(2:4), w(2:4)] = RK4(f, h, t(1), w(1));
    NFLAG = 1; % indicates computation from RK4
    i = 5; % CHANGED from 4 to 5 since indices start at 1, not 0
    t_temp = t(4) + h; % CHANGED from 3 to 4 since indices start at 1, not 0
    
    %disp("Step 4");
    while FLAG == 1
        %disp("Step 5");
        if (h == 0 || isnan(h))
            disp('WARNING: h is 0 or NaN!');
            disp(['h: ', num2str(h)]);
            break
        end
  
        WP = w(i-1) + (h/24) * (55*f(t(i-1), w(i-1)) - 59*f(t(i-2), w(i-2)) + 37*f(t(i-3), w(i-3)) - 9*f(t(i-4), w(i-4))); % predict w_i
        WC = w(i-1) + (h/24) * (9*f(t_temp, WP) + 19*f(t(i-1),w(i-1)) - 5*f(t(i-2), w(i-2)) + f(t(i-3), w(i-3))); % correct w_i
        sigma = 19*abs(WC-WP)/(270*h);
        %disp(['h: ', num2str(h)]); 
        %disp(['sigma: ', num2str(sigma)]);
        
        %disp("Step 6");
        if sigma <= TOL
            % do steps 7-16 (result accepted)
            %disp("Step 7");
            w(i) = WC;
            t(i) = t_temp;
            
            %disp("Step 8");
            if NFLAG == 1
                for j=(i-3):i
                    % OUTPUT (j, t_j, w_j, h)
                    disp([j-1 t(j) w(j) h]);
                    % previous results also accepted
                end
            else
                % OUTPUT (i, t_i, w_i, h)
                disp([i-1 t(i) w(i) h]);
                % previous results already accepted
            end
            
            %disp("Step 9");
            if LAST == 1
                FLAG = 0; % next step is 20
            else
                %disp("Step 10");
                i = i+1;
                NFLAG = 0;
                
                %disp("Step 11");
                if (sigma <= 0.1 * TOL || t(i-1) + h > b)
                    % increase h if it is more accurate than required, or
                    % decrease h to include b as a mesh point
                    %disp("Step 12");
                    q = (TOL / (2 * sigma)) ^ (1/4);
                    %disp(['h: ', num2str(h)]);
                    
                    %disp("Step 13");
                    if q > 4
                        h = 4*h;
                    else
                        h = q*h;
                    end
                    %disp(['h: ', num2str(h)]);
                    
                    %disp("Step 14");
                    if h > hmax
                        h = hmax;
                    end
                    %disp(['h: ', num2str(h)]);
                    
                    %disp("Step 15");
                    %disp([t(i-1), t(i-1) + 4*h, (b - t(i-1))/4]);
                    %disp(['h: ', num2str(h)]);
                    if (t(i-1) + 4*h > b)
                        %disp('TRUE');
                        h = (b - t(i-1))/4;
                        LAST = 1;
                    end
                    
                    
                    %disp("Step 16");
                    [t(i:i+2), w(i:i+2)] = RK4(f, h, t(i-1), w(i-1));
                    NFLAG = 1;
                    i = i+3;
                    %disp(['h: ', num2str(h)]);
                end
            end
            %disp(['h: ', num2str(h)]);
            
        else
            % do steps 17-19 (result rejected)
            %disp("Step 17");
            q = (TOL / (2 * sigma)) ^ (1/4);
            %disp(['h: ', num2str(h)]);
            
            %disp("Step 18");
            if q < 0.1
                h = 0.1*h;
            else
                h = q*h;
            end
            %disp(['h: ', num2str(h)]);
            
            %disp("Step 19");
            if h < hmin
                FLAG = 0;
                %disp("hmin exceeded");
            else
                if NFLAG == 1
                    i = i-3; % previous results also rejected
                end
                [t(i:i+2), w(i:i+2)] = RK4(f, h, t(i-1), w(i-1));
                i = i+3;
                NFLAG = 1; % end step 6
            end
        end
        
        %disp(['h: ', num2str(h)]);
        %disp("Step 20");
        t_temp = t(i-1) + h;
    end
end

% Helper function for Adams Variable SS.PC
function [y, z] = RK4(f, h, x, v)
    for j=2:4
        K1 = h * f( x(j-1), v(j-1) );
        K2 = h * f( x(j-1) + h/2, v(j-1) + K1/2 );
        K3 = h * f( x(j-1) + h/2, v(j-1) + K2/2 );
        K4 = h * f( x(j-1) + h, v(j-1) + K3);
        v(j) = v(j-1) + (K1 + 2*K2 + 2*K3 + K4) / 6;
        x(j) = x(1) + (j-1)*h;
    end
    y = x(2:4);
    z = v(2:4);
end