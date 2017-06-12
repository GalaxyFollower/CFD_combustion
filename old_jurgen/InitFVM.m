%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lambda=1 und D=1!!!!!!!

%%%%%%%%%%%%%% GEOMETRY%%%%%%%%%%%%%%

%    |---
%    |   ---
%    |      ----
%    |          |
% h1 |----------|  h2  <- symmetry axis
%    |          |
%    |      ----
%    |   ---
%    |---
%
%    |<--  l -->|


% Defin dimension of the trapezoidal domain
% h2 <= h1 !


shape = 'quadratic';  % 'linear' or 'quadratic'

h1 = 4;
hm = 4;            % only necessary for quatratic option 
h2 = 2;
l = 4;

% Number of degrees of freedom (number of nodes per length)

dimX = 12;
dimY = 10;


% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 2;
Tinf = 90;

% Boundary conditions (Only Dirichlet applied in Session 03) 

boundary.south = 'Neumann';
boundary.north = 'Robin';
boundary.east = 'Robin';
boundary.west = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 50;
TD.south= 50;
TD.west= 100;
TD.east=50;



% Shape of the Cooling Fin

% h2 <= h1 !

switch shape
    
    case 'linear'

        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'

        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;

        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;

end

%%%%%%%%%%%SCHEME%%%%%%%%%%%%%
Scheme='Steady'; %Steady or Unsteady Forward Euler or Unsteady Theta or Unsteady RK4
theta=0.5; %for theta between [0.5,1] scheme unconditionally stable
         %for theta between [0,0.5] scheme conditionally stable
T_u=ones(dimY,dimX)*Tinf;
switch Scheme
    case 'Steady'
        delta_t=0; %nicht wichtig
    case 'Unsteady Forward Euler'
        delta_t= 1/2*(l/dimX)^2*(h2/dimY)^2/((l/dimX)^2+(h2/dimY)^2)/4;%equation 4.39 with safety factor 1/4 (1/3 gives instability)
    case 'Unsteady Theta'
        if theta<0.5
            delta_t=1/(2*(1-2*theta))*(l/dimX)^2*(h2/dimY)^2/((l/dimX)^2+(h2/dimY)^2)/3;%equation 4.39 with safety factor 1/3
        else
            delta_t=0.01;
        end
    case 'Unsteady RK4'
        delta_t=1/2*(l/dimX)^2*(h2/dimY)^2/((l/dimX)^2+(h2/dimY)^2)/3;%equation 4.39 with safety factor 1/10 for sides!!!!
end

        