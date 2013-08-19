function [Uh, Vh, U, V] = FluidSolve(U, V, Fx, Fy, rho, mu, dx, dy, dt, Nx, Ny, Lx, Ly, MatDx, MatDy, MatLap, IndX, IndY)
% 
% This function solves the incompressible Navier-Stokes (NS) equations
%      rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + Fx
%      rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + Fy
%      u_x + v_y = 0.
%
% We are using Peskin's two-step algorithm where the advection terms
% is expressed in skew symmetric form.
%
% INPUTS:  
%	       U, V           The x and y components of the fluid velocity.  
%                         These should be input as matrices. 
%          Fx,Fy          The x and y components of the forcing on the
%                         grid.  These are matrices the same size as U and
%                         V.
%          rho            The fluid density, which is a scalar.
%          mu             The diffusion coefficient, which is a scalar.
%          dx, dy, dt     The spatial step-size. The value for a FULL time-step.
%          Nx, Ny         The number of grid points along each side of the
%                         domain.
%          Lx, Ly         The length of the domain.
%          MatDx, MatDy   Centered difference approximation to the first derivative 
%                         in the x and y direction with periodic BC.
%          MatLap         Centered 5-point laplacian stencil
%          IndX, IndY     Matrix containing indices in x- and y-direction.
%                         Used as the wave number in FFT solution.
%
% OUTPUTS: Uh, Vh, U, V   The updated x and y components of the fluid
%                         velocity at a full and half time-step.  
%                         These will be in matrix format.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

%
% Evolve the Fluid half a time-step
%

% Construct FFT Operator
A_hat = 1 + 2*mu*dt/rho*( (sin(pi*IndX/Nx)/dx).^2 + (sin(pi*IndY/Ny)/dy).^2 );

% Vectorize and construct Diagonal Matrix
UVec = U(:);
VVec = V(:);
UMat = spdiags(UVec, [0], length(UVec), length(UVec));
VMat = spdiags(VVec, [0], length(VVec), length(VVec));

% Construct right hand side in linear system
rhs_u = .5*dt/rho*( - .5*rho*(UMat*MatDx*UVec + VMat*MatDy*UVec) ...
                    - .5*rho*(MatDx*(UVec.^2) + MatDy*(VVec.*UVec)) ...
                    + Fx(:) ...
        ) + UVec;

rhs_v = .5*dt/rho*( - .5*rho*(UMat*MatDx*VVec + VMat*MatDy*VVec) ...
                    - .5*rho*(MatDx*(UVec.*VVec) + MatDy*(VVec.^2)) ...
                    + Fy(:) ...
        ) + VVec;

rhs_u = reshape(rhs_u,Ny,Nx);
rhs_v = reshape(rhs_v,Ny,Nx);
    
% Perform FFT
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  

% Calculate Fluid Pressure
p_hat = -(1.0i/dx*sin(2*pi*IndX/Nx).*rhs_u_hat + 1.0i/dy*sin(2*pi*IndY/Ny).*rhs_v_hat)./...
            ( 0.5*dt/rho*((sin(2*pi*IndX/Nx)/dx).^2 + (sin(2*pi*IndY/Ny)/dy).^2 ));
        
% Zero out modes.
p_hat(1,1) = 0;
p_hat(1,Nx/2+1) = 0;
p_hat(Ny/2+1,Nx/2+1) = 0;  
p_hat(Ny/2+1,1) = 0;

% Calculate Fluid Velocity
u_hat = (rhs_u_hat - 0.5*i*dt/(rho*dx)*sin(2*pi*IndX/Nx).*p_hat)./A_hat;
v_hat = (rhs_v_hat - 0.5*i*dt/(rho*dy)*sin(2*pi*IndY/Ny).*p_hat)./A_hat;

% 2d inverse of Fourier Coefficients
Uh = real(ifft2(u_hat));
Vh = real(ifft2(v_hat));

%
% Evolve the Fluid a full time-step
%

% FFT Operator is the same as in the first step.

% Vectorize and construct Diagonal Matrix
UhVec = Uh(:);
VhVec = Vh(:);
UhMat = spdiags(UhVec, [0], length(UhVec), length(UhVec));
VhMat = spdiags(VhVec, [0], length(VhVec), length(VhVec));

% Construct right hand side in linear system
rhs_u = dt/rho*( - .5*rho*(UhMat*MatDx*UhVec + VhMat*MatDy*UhVec) ...
                    - .5*rho*(MatDx*(UhVec.^2) + MatDy*(VhVec.*UhVec)) ...
                    + Fx(:) + mu/2*MatLap*UVec ...
        ) + UVec;

rhs_v = dt/rho*( - .5*rho*(UhMat*MatDx*VhVec + VhMat*MatDy*VhVec) ...
                    - .5*rho*(MatDx*(UhVec.*VhVec) + MatDy*(VhVec.^2)) ...
                    + Fy(:) + mu/2*MatLap*VVec ...
        ) + VVec;

rhs_u = reshape(rhs_u,Ny,Nx);
rhs_v = reshape(rhs_v,Ny,Nx);
    
% Perform FFT
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  

% Calculate Fluid Pressure
p_hat = -(1.0i/dx*sin(2*pi*IndX/Nx).*rhs_u_hat + 1.0i/dy*sin(2*pi*IndY/Ny).*rhs_v_hat)./...
            ( dt/rho*((sin(2*pi*IndX/Nx)/dx).^2 +(sin(2*pi*IndY/Ny)/dy).^2 ));
        
% Zero out modes.
p_hat(1,1) = 0;
p_hat(1,Nx/2+1) = 0;
p_hat(Ny/2+1,Nx/2+1) = 0;  
p_hat(Ny/2+1,1) = 0;

% Calculate Fluid Velocity
u_hat = (rhs_u_hat - 1.0i*dt/(rho*dx)*sin(2*pi*IndX/Nx).*p_hat)./A_hat;
v_hat = (rhs_v_hat - 1.0i*dt/(rho*dy)*sin(2*pi*IndY/Ny).*p_hat)./A_hat;

% 2d inverse of Fourier Coefficients
U = real(ifft2(u_hat));
V = real(ifft2(v_hat));

