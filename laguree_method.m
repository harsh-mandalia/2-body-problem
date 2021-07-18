% Laguerre?s Iteration Method Algorithm
clearvars;
iterations=zeros(100,315);  
E_array=zeros(100,315);
n=10;
for e = 0.01:0.01:0.99                  % Range of eccentricity 
    j=int16(e*100)+1;                   % for matrix
    for M = 0.01:0.01:3.14              % Range of Mean Anomaly
        k=int16(M*100)+1;
        E=M;                            % Initial Guess
        f_E = @(E) E-M-e*sin(E);        % Keplers Eqn function
        fdash = @(E) 1-e*cos(E);        % First derivative of function
        fddash = @(E) e*sin(E);
        i=0;
        % Laguerre's Formula
        while (true)
            g=(fdash(E))/f_E(E);
            h=g^2-((fddash(E))/f_E(E));
            a=n/(g+sign(g)*sqrt((n-1)*(n*h-g^2)));
            E=E-a;
            a/E;
            if(abs(a/E)<1e-4)
                break
            end
            i=i+1;   
            E_array(j,k)=E;
        end
        iterations(j,k)=i;
    end
end
surf(iterations)        % 3D plot of M,e,iterations
% surf(E_array)         % uncomment this line to see the 3d plot of M,e,E
shading interp
xlabel('M*100')
ylabel('e*100')
zlabel('number of iterations')
title("Number of iterations required for convergence using Laguerre's method")