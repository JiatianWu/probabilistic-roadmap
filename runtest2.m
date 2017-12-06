function[numofmoves, caught] = runtest(mapfile, armstart, armgoal, planner_id)

LINKLENGTH_CELLS=10;
envmap = load(mapfile);

close all;

%draw the environment
image(envmap'*255);
hold on;

t0 = clock;
%armplan should be a matrix of D by N 
%where D is the number of DOFs in the arm (length of armstart) and
%N is the number of steps in the plan 
armplan = armplanner(envmap, armstart, armgoal, planner_id);

time = etime(clock, t0);

finalplan = armplan;
%interpolation
interp = 1;
if (interp == 1)
    finalplan = [];
    for j = 1:size(armplan, 1) - 1   %waypoint
        distance = 0;
        for i = 1:size(armplan, 2)  %DOFs
            if (distance < abs(angdiff(armplan(j+1, i), armplan(j, i))))
                distance = abs(angdiff(armplan(j+1, i), armplan(j, i)));
            end
        end
        steps = fix(distance / (pi / 20));
        if (steps < 2)
           steps = steps;
        end
        interplan = [];
        for i = 1:size(armplan, 2)  %DOFs
            for k = 1:steps
                interplan(k,i) = armplan(j,i) + ((k-1)/(steps-1) * (angdiff(armplan(j, i), armplan(j+1, i))));
            end
        end
        if size(interplan) ~= 0
    %     finalplan = horzcat(finalplan, interplan);
          finalplan = [finalplan; interplan];
        end
    end
end

    
fprintf(1, 'plan of time %f \n', time);
fprintf(1, 'plan of length %d was found\n', size(armplan,1));

%draw the plan

midx = size(envmap,2)/2;
x = zeros(length(armstart)+1,1);
x(1) = midx;
y = zeros(length(armstart)+1,1);
for i = 1:size(finalplan)
    for j = 1:length(armstart)
        x(j+1) = x(j) + LINKLENGTH_CELLS*cos(finalplan(i,j));
        y(j+1) = y(j) + LINKLENGTH_CELLS*sin(finalplan(i,j));
    end;
    plot(x,y, 'c-');
    pause(0.2);
end;

%armplan