function[] = main(map)

angles = [pi/2 pi/4 pi/2 pi/4 pi/2;
          pi/8 3*pi/4 pi 0.9*pi 1.5*pi;
          pi/2 pi/2 pi/2 pi/2 0;
          pi/3 pi/3 0.9*pi pi/2 pi/2];
 startQ = angles(1,:);  
 goalQ = angles(2,:);
 runtest('map1.txt', startQ, goalQ, 3)
      
% pathSum = 0;
% timeSum = 0;
% pathMean = 0;
% timeMean = 0;
% if map == 1
%     for i = 1:4
%         startQ = angles(i,:);
%         for j = 1:4
%             if j ~= i
%                 goalQ = angles(j,:);
%                 [path, time] = runtest('map1.txt', startQ, goalQ, 1);
%                 pathSum = pathSum + path;
%                 timeSum = timeSum + time;
%             end;
%         end;
%     end;
%     pathMean = pathSum / 12;
%     timeMean = timeSum / 12;
%     fprintf(2, 'path %f time %f \n', pathMean, timeMean);
% elseif map == 2
%     for i = 1:4
%         startQ = angles(i,:);
%         for j = 1:4
%             if j ~= i
%                 goalQ = angles(j,:);
%                 [path, time] = runtest('map2.txt', startQ, goalQ, 1);
%                 pathSum = pathSum + path;
%                 timeSum = timeSum + time;
%             end;
%         end;
%     end;   
%     pathMean = pathSum / 12;
%     timeMean = timeSum / 12;
% end    