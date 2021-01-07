function nbody_conv_test()
% Tests the nbody function and does a 3-level convergence test

    % set params for nbody
    tmax = 100;
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
    
    levels = [6,7,8,9];
    [t6, r6, v6] = nbody(tmax, levels(1), [m1; m2], r0, v0);
    [t7, r7, v7] = nbody(tmax, levels(2), [m1; m2], r0, v0);
    [t8, r8, v8] = nbody(tmax, levels(3), [m1; m2], r0, v0);
    [t9, r9, v9] = nbody(tmax, levels(4), [m1; m2], r0, v0);
    
    % get x coordinate from particle 1
    r6 = squeeze(r6(1,1,:));
    r7 = squeeze(r7(1,1,:));
    r8 = squeeze(r8(1,1,:));
    r9 = squeeze(r9(1,1,:));
    
    % make plot to show deviations
    clf;
    figure(1);
    hold on;
    plot(t6, r6, 'r-.o');
    plot(t7, r7, 'g-.+');
    plot(t8, r8, 'b-.*');
    xlabel("time")
    ylabel("x coordinate")
    legend(['level=6'; 'level=7'; 'level=8']);
    
    % downsample data generated with level 7 and 8
    r7 = r7(1:2:end);
    r8 = r8(1:4:end);
    r9 = r9(1:8:end);
    
    dr67 = r6 - r7;
    dr78 = r7 - r8;
    dr89 = r8 - r9;
    
    % make plot to show differences between different levels
    figure(2);
    hold on;
    plot(t6, dr67, 'r-.o');
    plot(t6, dr78, 'g-.+');
    plot(t6, dr89, 'b-.+');
    xlabel("time");
    legend(["r6 - r7", "r7 - r8", "r8 - r9"], 'Location', 'northwest');
    
    % scale dr78 by 4
    dr78 = 4 * dr78;
    dr89 = 16 * dr89;
    figure(3);
    hold on;
    plot(t6, dr67, 'r-.o');
    plot(t6, dr78, 'g-.+');
    plot(t6, dr89, 'b-.+');
    xlabel("time");
    legend(["r6 - r7", "4*(r7 - r8)", "16*(r8 - r9)"], 'Location', 'northwest');
    
end
    
