% Newton's Iteration Method Algorithm
clearvars;
iterations=zeros(100,315); 
E_array=zeros(100,315);
for e = 0.01:0.01:0.99                  % Range of eccentricity 
    j=int16(e*100)+1;                   % for matrix
    for M = 0.01:0.01:3.14              % Range of Mean Anomaly
        k=int16(M*100)+1;
        E=M;                            % Initial Guess
        f_E = @(E) E-M-e*sin(E);        % Keplers Eqn function
        fdash = @(E) 1-e*cos(E);        % First derivative of function
        i=0;
        % Newton's Formula
        while (true)    
            a=f_E(E)/fdash(E);
            E = E - a;   
            if(abs(a/E)<1e-4)
                break
            end
            i=i+1;                                         
        end
        E_array(j,k)=E;
        iterations(j,k)=i;
    end
end
surf(iterations)        % 3D plot of M,e,iterations
% surf(E_array)         % uncomment this line to see the 3d plot of M,e,E
shading interp
xlabel('M*100')
ylabel('e*100')
zlabel('number of iterations')
title("Number of iterations required for convergence using Newton's method")