function [t, r, v] = nbody(tmax, level, mass, r0, v0)
% Solves nbody problem using second order finite difference approximations
%
% Input arguments
%
%     tmax:   (real scalar) Final solution time
%     level:  (integer scalar) Discretization level
%     mass:   (N x 1 array) Particle masses
%     r0:     (N x 3 array) Intial particle position vectors
%     v0:     (N x 3 array) Initial particle velocity vectors
%  
% Return value
% 
%      t:     (real vector) Vector of length nt = 2^level + 1 containing
%             discrete times (time mesh)
%      r:     (N x 3 x nt) positions of all particles at all time steps
%      v:     (N x 3 x nt) velocities of all particles at all time steps

   % trace controls "tracing" output.  Set 0 to disable, non-0 to enable.
   trace = 1;
   % tracefreq controls frequency of tracing output in main time step loop.
   tracefreq = 100;

   if trace
      fprintf('In nbody: Argument dump follows\n');
      tmax, level, mass, r0, v0
   end
   
   s = size(mass);
   N = s(1);

   % Define number of time steps and create t, theta and omega arrays of
   % appropriate size for efficiency.
   nt = 2^level + 1;
   t = linspace(0.0, tmax, nt);
   r = zeros(N, 3, nt);
   v = zeros(N, 3, nt);
   
   % Determine discrete time step from t array.
   deltat = t(2) - t(1);
   
   % Initilize the initial values of r and v
   r(:,:,1) = r0;
   v(:,:,1) = v0;
   
   % Initilize the second value for r
   r(:,:,2) = r0 + v0 * deltat + 0.5 * deltat^2 * nbodyaccn(mass, r0);
   
  if trace
      fprintf('deltat=%g, r1, r2\n', deltat);
      r1 = r(:,:,1)
      r2 = r(:,:,2)
  end
  
   % Evolve the nbody problem to the final time using the discrete equations
   % of motion.  Also compute an estimate of the velocity vector at
   % each time step.
   for n = 2 : nt-1
       % This generates tracing output every 'tracefreq' steps.
       if rem(n, tracefreq) == 0
           fprintf('nbody: Step %d of %d\n', n, nt);
       end
       
       r(:,:,n+1) = deltat^2 * nbodyaccn(mass, r(:,:,n)) + 2 * r(:,:,n) - r(:,:,n-1);
       
       v(:,:,n) = (r(:,:,n+1) - r(:,:,n-1)) / 2 * deltat;
   end
   % Use linear extrapolation to determine value of v at t = nt
   v(:,:,nt) = 2 * v(:,:,nt-1) - v(:,:,nt-2);
end
   