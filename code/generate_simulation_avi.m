function generate_simulation_avi(avifilename, timesteps, r)
% Produces an animation of the galaxy dynamics/collisions
%
% Input arguments
%
%     avifilename: (string) filename of the avi video
%     timesteps:   (1 X nt array) Array containing all timesteps
%     r:           (N x 3 x nt array) Vector containing positions of all
%                                     particles at all timesteps

% Setting this parameter to a "largish" value, say 0.1
% (seconds), will produce a slow-motion effect.
% 
% Set it to 0 for maximum animation speed.
pausesecs = 0.0;

% Plot attributes defining the appearance of the stars.
star_size = 1;
star1_color = 'm';
star2_color = 'c';
star_marker = 'o';

% Compute size of plot
x = r(:,1,1:end/2);
pos_x = mean(x(x>0));
neg_x = mean(x(x<0));

y = r(:,2,1:end/2);
pos_y = mean(y(y>0));
neg_y = mean(y(y<0));

lmax = max(pos_x, pos_y);
lmin = min(neg_x, neg_y);
mult = 5;

% Presumed AVI playback rate in frames per second.
aviframerate = 25;

% init video write object
aviobj = VideoWriter(avifilename, 'MPEG-4');
open(aviobj);

% number of particles
s = size(r);
num_stars = (s(1)-2)/2;

% number of timesteps
s_t = size(timesteps);
num_timesteps = s_t(2);

for t = 1:num_timesteps 
      fprintf("timestep # %d\n", t);
      % Clear figure
      clf;

      % Don't erase figure after each plot command.
      hold on;
      
      % set background color to black
      set(gca,'Color','k')
      % remove axes 
      set(gca,'xtick',[])
      set(gca,'ytick',[]);
      set(gca,'Units','pixels','Position',[0 0 565 565])

      % Define plotting area, using a 1:1 aspect ratio for the 
      % plotted region, boxed axes and a 15%-width "border" around 
      % the unit square.
      axis square;
      xlim([mult*lmin, mult*lmax]);
      ylim([mult*lmin, mult*lmax]);
      
      % draw the particles of the first core
      plot(r(3:num_stars+2, 1, t), r(3:num_stars+2, 2, t), 'Marker', star_marker, 'MarkerSize', star_size, ...
           'MarkerEdgeColor', star1_color, 'MarkerFaceColor', star1_color, 'LineStyle', 'none');
      
      % draw the particles of the second core
      % draw the particle
      plot(r(3+num_stars:end, 1, t), r(3+num_stars:end, 2, t), 'Marker', star_marker, 'MarkerSize', star_size, ...
           'MarkerEdgeColor', star2_color, 'MarkerFaceColor', star2_color, 'LineStyle', 'none');

      % Force update of figure window.
      drawnow;

      % Record video frame and record multiple copies of the first frame. 
      % Here we record 5 seconds worth which will allow the viewer a bit
      % of time to process the initial scene before the animation starts.
      if t == 0
          framecount = 5 * aviframerate;
      else
          framecount = 1;
      end
      for iframe = 1 : framecount
          writeVideo(aviobj, getframe(gcf));
      end

      % Pause execution to control interactive visualization speed.
      pause(pausesecs);
end

close(aviobj);
fprintf('Created video file: %s\n', avifilename);
end
