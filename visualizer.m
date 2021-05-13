function [] = visualizer()

    close all;
    figure; hold on;

    % change these to test different maps and start/goal configs
    mapChoice = 3;
    startGoal = 3;
    % change line 159 to adjust the default view in plot3()
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % draw the space with obsatcles
    gridsize = 100;
    gridmeters = 5;
    x = linspace(0,gridmeters,gridsize);
    y = x;
    zb = zeros(gridsize,gridsize);
    zt = gridmeters*ones(gridsize,gridsize);
    
    % obstacles (yrange, xrange)
    % map 1
    if mapChoice == 1
        zb(10:20,10:25) = gridmeters;
        zb(38:50,79:83) = gridmeters;

        zt(25:50,50:60) = 0.25;
        zt(43:70,10:16) = 0.25;
        zt(70:90,30:55) = 0.25;
    end
    
    % map 2
    if mapChoice == 2
        zt(10:20,10:25) = 0.25;
        zt(38:50,79:83) = 0.25;
        
        zb(25:50,50:60) = gridmeters;
        zb(43:70,10:16) = gridmeters;
        zb(70:90,30:55) = gridmeters;
    end
    
    % map 3
    if mapChoice == 3
        zb(1:50,10:30) = gridmeters;
        zt(51:60,10:30) = 0.25;
        zb(61:85,10:30) = gridmeters;
    end
       
    % start/goal 1
    if startGoal == 1
        armstart = [0,0,0,0,0.25,3,0];
        armgoal = [pi/2,pi/2,-pi/2,pi/2,4.5,4.5,pi];
    end
    
    % start/goal 2
    if startGoal == 2
        armstart = [0,pi/2,0,0,0.25,0.25,1];
        armgoal = [pi/2,pi/2,pi/2,-pi/2,2,4.75,pi/4];
    end
    
    % start/goal 3
    if startGoal == 3
        armstart = [0,pi/2,0,0,0.05,0.05,0];
        armgoal = [pi/2,pi/2,-pi/2,pi/2,4.5,0.1,0];
    end
        

    surf(x,y,zb);
    surf(x,y,zt);
    xlim([0,gridmeters]); ylim([0,gridmeters]); zlim([0,gridmeters]);

    % create a 3D array from zb, zt to pass to planner
    map = zeros(gridsize,gridsize,gridsize);
    for heightIndex = 1:gridsize
        height = heightIndex*gridmeters/gridsize;
        for row = 1:gridsize
            for col = 1:gridsize
               if (zb(row,col) >= height || zt(row,col) <= height)
                   map(row,col,heightIndex) = 1;
               end
            end
        end
    end
    % armstart = [0,0,0,0,0.25,3,0];
    % armgoal = [pi/2,pi/2,-pi/2,pi/2,4.5,4.5,pi];
    % map(:,:,50)
    planner_id = 1;
%     map = zeros(3,4,2);
%     map(1,:,1) = [1 2 3 4];
%     map(2,:,1) = [5,6,7,8];
%     map(3,:,1) = [9,10,11,12];
%     map(1,:,2) = [13,14,15,16];
%     map(2,:,2) = [17,18,19,20];
%     map(3,:,2) = [21,22,23,24];
    
    % pass the env to the planner
    % envmap --> z, 2D to 3D array
    % [armplan, armplanlength, nodesgenerated, plantime] = armplanner(envmap, armstart, armgoal, planner_id); 
    [armplan, armplanlength, nodesgenerated, plantime] = armplanner(map, armstart, armgoal, planner_id); 
    % run the robot
    l0 = 0.1;
    l1 = 0.35;
    l2 = 0.35;
    l3 = 0.125;

%     n = 100;
%     th0 = linspace(0, pi/2, n);
%     th1 = linspace(0, pi/2, n);
%     th2 = linspace(0, -pi/2, n);
%     th3 = linspace(0, pi/2, n);
% 
%     xpos = linspace(0.25, 4.5, n);
%     ypos = linspace(3, 4.5, n); % sin(xpos);
%     thpos = linspace(0,pi,n);

    n = size(armplan,1);
    th0 = armplan(:,1);
    th1 = armplan(:,2);
    th2 = armplan(:,3);
    th3 = armplan(:,4);
    
    xpos = armplan(:,5);
    ypos = armplan(:,6);
    thpos = armplan(:,7);

    for i = 1:n

        r0 = 0; 
        r1 = 0;
        r2 = l1*cos(th1(i));
        r3 = r2 + l2*cos(th1(i)+th2(i));
        r4 = r3 + l3*cos(th1(i)+th2(i)+th3(i));

        x0 = 0;
        x1 = 0;
        x2 = (r2-r1)*cos(th0(i)+thpos(i));
        x3 = (r3-r1)*cos(th0(i)+thpos(i));
        x4 = (r4-r1)*cos(th0(i)+thpos(i));

        y0 = 0;
        y1 = 0;
        y2 = (r2-r1)*sin(th0(i)+thpos(i));
        y3 = (r3-r1)*sin(th0(i)+thpos(i));
        y4 = (r4-r1)*sin(th0(i)+thpos(i));

        z0 = 0;
        z1 = z0 + l0;
        z2 = z1 + l1*sin(th1(i));
        z3 = z2 + l2*sin(th1(i) + th2(i));
        z4 = z3 + l3*sin(th1(i) + th2(i) + th3(i));

        x = [x0 x1 x2 x3 x4] +xpos(i);
        y = [y0 y1 y2 y3 y4] +ypos(i);
        z = [z0 z1 z2 z3 z4];

        plot3([x0,x1]+xpos(i),[y0,y1]+ypos(i),[z0,z1],'-r');
        plot3([x1,x2]+xpos(i),[y1,y2]+ypos(i),[z1,z2],'-g');
        plot3([x2,x3]+xpos(i),[y2,y3]+ypos(i),[z2,z3],'-b');
        plot3([x3,x4]+xpos(i),[y3,y4]+ypos(i),[z3,z4],'-c');
        view([-50, 10]); % adjust this if desired
        pause(0.1);

    end

end