function [Uh, Vh, U, V] = FluidSolve(U, V, Fx, Fy, rho, mu, dx, dy, dt, Nx, Ny, Lx, Ly, IndX, IndY)
% 
% This function solves the incompressible Navier-Stokes (NS) equations
%      rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + Fx
%      rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + Fy
%      u_x + v_y = 0.
%
% We are using Peskin's two-step algorithm where the advection terms
% is expressed in skew symmetric form.
%
% INPUTS:  U, V           The x and y components of the fluid velocity.  
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

% Compute first and second derivatives (centred).
Ux = Centered_Dx(U,dx);
Uy = Centered_Dy(U,dy);
Uxx = Centered_Dxx(U,dx);
Uyy = Centered_Dyy(U,dy);
Vx = Centered_Dx(V,dx);
Vy = Centered_Dy(V,dy);
Vxx = Centered_Dxx(V,dx);
Vyy = Centered_Dyy(V,dy);

% Computed derivatives of products U^2, V^2, and U*V.
USq = U.^2;
VSq = V.^2;
UV = U.*V;
USq_x = Centered_Dx(USq,dx);
VSq_y = Centered_Dy(VSq,dy);
UV_x = Centered_Dx(UV,dx);
UV_y = Centered_Dy(UV,dy);

% Construct right hand side in linear system
rhs_u = .5*dt/rho*( - .5*rho*(U.*Ux + V.*Uy) ...
                    - .5*rho*(USq_x + UV_y) ...
                    + Fx ...
        ) + U;

rhs_v = .5*dt/rho*( - .5*rho*(U.*Vx + V.*Vy) ...
                    - .5*rho*(UV_x + VSq_y) ...
                    + Fy ...
        ) + V;

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

% Compute first derivatives (centred) at half step.
Uhx = Centered_Dx(Uh,dx);
Uhy = Centered_Dy(Uh,dy);
Vhx = Centered_Dx(Vh,dx);
Vhy = Centered_Dy(Vh,dy);

% Computed derivatives of products U^2, V^2, and U*V at half step.
UhSq = Uh.^2;
VhSq = Vh.^2;
UhVh = Uh.*Vh;
UhSq_x = Centered_Dx(UhSq,dx);
VhSq_y = Centered_Dy(VhSq,dy);
UhVh_x = Centered_Dx(UhVh,dx);
UhVh_y = Centered_Dy(UhVh,dy);

% Construct right hand side in linear system
rhs_u = dt/rho*( - .5*rho*(Uh.*Uhx + Vh.*Uhy) ...
                    - .5*rho*(UhSq_x + UhVh_y) ...
                    + Fx + mu/2*(Uxx+Uyy) ...
        ) + U;

rhs_v = dt/rho*( - .5*rho*(Uh.*Vhx + Vh.*Vhy) ...
                    - .5*rho*(UhVh_x + VhSq_y) ...
                    + Fy + mu/2*(Vxx+Vyy) ...
        ) + V;

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

