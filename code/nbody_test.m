function nbody_test()
% Tests the nbody function

    % set params for nbody
    tmax = 40.0;
    level = 8;
    m1 = 2;
    m2 = 5;
    % set params according to the suggested test case
    m = m1 + m2;
    r = 4;
    r1 = m2/m*r;
    r2 = m1/m*r;
    % initialize the r0 vector
    r0 = zeros(2, 3);
    r0(1,:) = [r1 0 0];
    r0(2,:) = [-r2 0 0];
    % initialize the v0 vector
    v0 = zeros(2, 3);
    v0(1,:) = [0 sqrt(m2*r1)/r 0];
    v0(2,:) = [0 -sqrt(m1*r2)/r 0];
    
    % call nbody
    [t, r, v] = nbody(tmax, level, [m1; m2], r0, v0);
    clf; hold on;
    %grid;
    plot(squeeze(r(1,1,:)), squeeze(r(1,2,:)), 'r');
    plot(squeeze(r(2,1,:)), squeeze(r(2,2,:)), 'g');
    
    % draw arrows to show starting point + initial velocity vector
    s = 10;
    p1 = [r1 0];
    p2 = [r1 r(1,2,s)];
    p3 = [-r2 0];
    p4 = [-r2 r(2,2,s)];
    dp1 = p2-p1;
    dp2 = p4-p3;
    quiver(p1(1),p1(2),dp1(1),dp1(2),0)
    quiver(p3(1),p3(2),dp2(1),dp2(2),0)
    
    % add labels + legend to plot
    xlabel("x coordinate")
    ylabel("y coordinate")
    l1 = sprintf('particle 1, mass = %d', m1);
    l2 = sprintf('particle 2, mass = %d', m2);
    legend([l1; l2]);
    
end
    
