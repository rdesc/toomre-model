function galaxy_sim()
% simulates the dynamics of two galaxies with orbiting stars using the toomre model
    fprintf('==============================================\n')
 
    %%%%%%%%%%%%%%%%%%% main params of the simulation %%%%%%%%%%%%%%%%%%%%%
    gen_avi = 1; % option to generate avi from simulation
    avifilename = 'toomre.mp4';
    
    plot_t0 = 0; % option to plot star positions at t0
    plot_traj = 0; % option to plot star trajectories for all time steps
 
    % initialize number of stars per core
    num_stars = 5000;
    % initialize the rotation direction of the stars (1 = CW, -1 = CCW)
    rot1 = 1;
    rot2 = -1;
    % initialize max and min size of stellar radii from core
    min_rad = 0.5;
    max_rad = 2;
    % option whether to evenly distribute stars in all 3 dimensions
    three_d = -1; % if -1 then just distribute stars in xy plane
 
    tmax = 10.0;
    level = 8;

    % initialize the masses
    m1 = 10;
    m2 = 4;
    m = m1 + m2;
    % initialize the r0 vector
    r = 5;
    r1 = m2/m*r;
    r2 = m1/m*r;
    r0 = zeros(2, 3);
    r0(1,:) = [5 -5 0];
    r0(2,:) = [0 0 0];
%     r0(1,:) = [r1 0 0];
%     r0(2,:) = [-r2 0 0];

    % initialize the v0 vector
    v0 = zeros(2, 3);
%     v0(1,:) = [0 sqrt(m2*r1)/r 0];
%     v0(2,:) = [0 -sqrt(m1*r2)/r 0];
    v0(1,:) = [0 1.0 0];
    v0(2,:) = [0 -1.0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%% end of main params %%%%%%%%%%%%%%%%%%%%%%%%%

    % randomly generate stars positions
    [s1_r0, s2_r0] = generate_stars(num_stars, min_rad, max_rad, three_d);
    
    % update star positions based on core positions
    s1_r0 = s1_r0 + r0(1,:);
    s2_r0 = s2_r0 + r0(2,:);
    
    % initialize star velocities such that they orbit cores
    s1_v0 = init_star_vel(m1, r0(1,:), s1_r0, rot1);
    s2_v0 = init_star_vel(m2, r0(2,:), s2_r0, rot2);
    
    % update star velocities based on core initial velocities
    s1_v0 = s1_v0 + v0(1,:);
    s2_v0 = s2_v0 + v0(2,:);
    
    if plot_t0
        % plot the initial positions of the two galaxies
        clf;
        hold on;
        grid;
        scatter(s1_r0(:,1), s1_r0(:,2))
        scatter(s2_r0(:,1), s2_r0(:,2))
        xlabel("x coordinate")
        ylabel("y coordinate")
        %legend(["stars from galaxy 1"; "stars from galaxy 2"], 'Location', 'northwest');
        title("Initial galaxy positions")

        % update positions based on initial velocities
        deltat = 0.1;
        s1_r1 = s1_r0 + s1_v0 * deltat;
        s2_r1 = s2_r0 + s2_v0 * deltat;

        % plot the initial velocity vectors of the star
        % draw arrows to show starting point + initial velocity vector
        for i = 1:num_stars
            p1 = [s1_r0(i,1) s1_r0(i,2)];
            p2 = [s1_r1(i,1) s1_r1(i,2)];
            p3 = [s2_r0(i,1) s2_r0(i,2)];
            p4 = [s2_r1(i,1) s2_r1(i,2)];
            dp1 = p2-p1;
            dp2 = p4-p3;
            quiver(p1(1),p1(2),dp1(1),dp1(2),0)
            quiver(p3(1),p3(2),dp2(1),dp2(2),0)
        end
    end
    
    % init mass, r0, v0 for nbody args
    mass = [m1; m2; zeros(2*num_stars, 1)];
    f_r0 = [r0; s1_r0; s2_r0];
    f_v0 = [v0; s1_v0; s2_v0];

    [t, r, v] = nbody(tmax, level, mass, f_r0, f_v0);
    
    % plot all trajectories
    if plot_traj
        figure(2);
        hold on;
        plot(squeeze(r(1,1,:)), squeeze(r(1,2,:)), 'b');
        plot(squeeze(r(2,1,:)), squeeze(r(2,2,:)), 'g');
        for i = 3:2*num_stars + 2
            plot(squeeze(r(i,1,:)), squeeze(r(i,2,:)), 'r');
        end
        xlabel("x coordinate")
        ylabel("y coordinate")
        dim = [0.2 0.5 0.3 0.3];
        str = sprintf('# stars/core = %d, tmax = %d', num_stars, tmax);
        t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
        t.BackgroundColor = 'white';
        title("Star trajectories at all timesteps")
    end
    
    % generate avi from simulation results
    if gen_avi
        generate_simulation_avi(avifilename, t, r);
    end
    fprintf('Done simulation!');

function v0 = init_star_vel(c_mass, c_r0, s_r0, rot)
% Initializes stellar velocities such that they orbit their core
%
% Input arguments
%
%     c_mass: (real scalar) Mass of the core
%     c_r0:   (1 X 3 array) Initial position of the core x,y,z
%     s_r0:   (N x 3 array) Initial positions of the stars
%     rot:    (-1 or 1) Direction of the orbit 1=CW, -1=CCW
%  
% Return value
% 
%      v0:     (N x 3 arrat) Intial velocity vector for the stars

   s = size(s_r0);
   N = s(1);
   rot = rot * pi/2;
   
   % expand core dimension to match dimension of stars
   c_r0 = repmat(c_r0, N, 1);
   
   rij = sqrt((c_r0(:,1) - s_r0(:,1)).^2 + ...
              (c_r0(:,2) - s_r0(:,2)).^2 + ...
              (c_r0(:,3) - s_r0(:,3)).^2);
          
   % transform separation vector by +/-90 degrees
   Rz90 = [cos(rot) -sin(rot)  0; ...
           sin(rot)  cos(rot)  0; ...
           0         0         1];
   
   sep_v = zeros(N, 3);
   diff = (c_r0 - s_r0);
   for i=1:N
       sep_v(i,:) = transpose(Rz90 * transpose(diff(i,:)));
   end

   v0 = sqrt(c_mass ./ rij) .* sep_v ./ rij;


function [r1, r2] = generate_stars(num_stars, min_rad, max_rad, three_d)
% Generate randomly distributed stars
%
% Input arguments
%
%     num_stars: (real scalar) Number of stars for a single core
%     min_rad:   (real scalar) Minimum radius from core
%     max_rad:   (real scalar) Maximum radius from core
%     three_d:   (-1 or 1) whether to distribute stars evenly in all 3
%                          dimensions, if -1 then just xy plane
%  
% Return values
% 
%      r1: (num_stars x 3 array) Initial position vecotr for stars
%                                around core 1
%      r2: (num_stars x 3 array) Initial position vector for stars
%                                around core 2

    % generate random radius values between min and max
    rad1 = (max_rad - min_rad).*rand(num_stars,1) + min_rad;
    rad2 = (max_rad - min_rad).*rand(num_stars,1) + min_rad;
    
    % generate random angle theta values
    min_ = -pi;
    max_ = pi;
    theta1 = (max_ - min_).*rand(num_stars,1) + min_;
    theta2 = (max_ - min_).*rand(num_stars,1) + min_;
    
    % generate random angle phi values
    phi1 = (max_ - min_).*rand(num_stars,1) + min_;
    phi2 = (max_ - min_).*rand(num_stars,1) + min_;
    
    % compute euclidian coordinates
    if three_d == 1
        r1 = [rad1 .* sin(theta1) .* cos(phi1),...
              rad1 .* sin(theta1) .* sin(phi1),...
              rad1 .* cos(theta1)];
        r2 = [rad2 .* sin(theta2) .* cos(phi2),...
              rad2 .* sin(theta2) .* sin(phi2),...
              rad2 .* cos(theta2)];
    else
        % compute euclidian coordinates (ignore z direction)
        r1 = [rad1 .* cos(theta1), rad1 .* sin(theta1), zeros(num_stars, 1)];
        r2 = [rad2 .* cos(theta2), rad2 .* sin(theta2), zeros(num_stars, 1)];
    end
    