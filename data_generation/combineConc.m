function [F_conc,Points_conc] = combineConc(V_steelcorr_surf,Points_corr_surf,Points_rand,...
            F_conc_outer,Points_conc_outer,elemlength,radius,rebarR,F_corr_surf,conc_spacing)
        
        % find top and bottom of corrosion surf
        indexTopCorr = find(Points_corr_surf(:,3) == elemlength);
        indexBotCorr = find(Points_corr_surf(:,3) == 0);
        
        % remove surface indexing of points
        indicator1 = zeros(size(Points_conc_outer,1),1);
        indicator1(Points_conc_outer(:,3) > (elemlength-10^(-10)),1) = 1;
        indicator1(Points_conc_outer(:,3) < 10^(-10),1) = 1;
        indicator2 = indicator1(F_conc_outer);
        F_conc_outer = F_conc_outer(mean(indicator2,2) ~= 1,:);
        
        % remove points inside corrosion surf
        indexInConc = logical(1 - checkInCircleVec(Points_conc_outer(:,1),Points_conc_outer(:,2),...
            radius,radius,rebarR));
        Points_conc_outer = Points_conc_outer(indexInConc,:);
        
        % shift surf of conc outer to match indexing of new points vec
        indexInConc = double(indexInConc);
        indexInConc(indexInConc == 1) = 1:size(Points_conc_outer,1);
        F_conc_outer = indexInConc(F_conc_outer);
        
        % find top and bottom of concrete outer surf
        indexTopConc = find(Points_conc_outer(:,3) > (elemlength-10^(-10)));
        indexBotConc = find(Points_conc_outer(:,3) < 10^(-10));
        Points_conc_outer(indexTopConc,3) = elemlength;
        Points_conc_outer(indexBotConc,3) = 0;
        
        % find outer circle of points
        [thetaTop,rhoTop] = cart2pol(Points_conc_outer(indexTopConc,1)-radius,Points_conc_outer(indexTopConc,2)-radius);
        [thetaBot,rhoBot] = cart2pol(Points_conc_outer(indexBotConc,1)-radius,Points_conc_outer(indexBotConc,2)-radius);
        indexTopOutConc = find(abs(rhoTop) > (radius-10^(-10)) & Points_conc_outer(indexTopConc,3) > (elemlength-10^(-10)));
        indexBotOutConc = find(abs(rhoBot) > (radius-10^(-10)) & Points_conc_outer(indexBotConc,3) < 10^(-10));
%         indexTopOutConc = indexTopConc(indexTopOutConc,1);
%         indexBotOutConc = indexBotConc(indexBotOutConc,1);
        
        % sort outer points of concrete
        [~,sortIndexTopOutConc] = sort(thetaTop(indexTopOutConc,1));
        [~,sortIndexBotOutConc] = sort(thetaBot(indexBotOutConc,1));
        indexTopOutConc = indexTopOutConc(sortIndexTopOutConc,:);
        indexBotOutConc = indexBotOutConc(sortIndexBotOutConc,:);
        
        % delaunay triangulate top and bottom
        Points_Top = [Points_corr_surf(indexTopCorr,1:2); Points_conc_outer(indexTopConc,1:2)];  
        C_top = [(1:(size(indexTopCorr,1)-1))' (2:size(indexTopCorr,1))';
            size(indexTopCorr,1) 1;
            indexTopOutConc(1:(end-1),1)+size(indexTopCorr,1) indexTopOutConc(2:end,1)+size(indexTopCorr,1);
            indexTopOutConc(end,1)+size(indexTopCorr,1) indexTopOutConc(1,1)+size(indexTopCorr,1)];
        DT_top = delaunayTriangulation(Points_Top,C_top);
        
        Points_Bot = [Points_corr_surf(indexBotCorr,1:2);Points_conc_outer(indexBotConc,1:2)];
        C_bot = [(1:(size(indexBotCorr,1)-1))' (2:size(indexBotCorr,1))';
            size(indexBotCorr,1) 1;
            indexBotOutConc(1:(end-1),1)+size(indexBotCorr,1) indexBotOutConc(2:end,1)+size(indexBotCorr,1);
            indexBotOutConc(end,1)+size(indexBotCorr,1) indexBotOutConc(1,1)+size(indexBotCorr,1)];
        DT_Bot = delaunayTriangulation(Points_Bot,C_bot);
        
        % remove middle points
        X = DT_top.Points(:,1);
        X_top = mean(X(DT_top.ConnectivityList),2);
        Y = DT_top.Points(:,2);
        Y_top = mean(Y(DT_top.ConnectivityList),2);
        
        inCircleTop = logical(1 - checkInCircleVec(X_top,Y_top,radius,radius,rebarR));
        F_top = DT_top.ConnectivityList(inCircleTop,:);
        
        X = DT_Bot.Points(:,1);
        X_bot = mean(X(DT_Bot.ConnectivityList),2);
        Y = DT_Bot.Points(:,2);
        Y_bot = mean(Y(DT_Bot.ConnectivityList),2);
        
        inCircleBot = logical(1 - checkInCircleVec(X_bot,Y_bot,radius,radius,rebarR));
        F_bot = DT_Bot.ConnectivityList(inCircleBot,:);
        
        % creat complete vector of points and surface
%         Points_conc = [V_steelcorr_surf; Points_conc_outer; Points_rand];
        Points_conc = [Points_corr_surf; Points_conc_outer; Points_rand];
        
        F_top_shift = zeros(size(F_top,1),3);
        F_top_shift(F_top(:,1) <= size(indexTopCorr,1),1) =...
            indexTopCorr(F_top(F_top(:,1) <= size(indexTopCorr,1),1),1);
        F_top_shift(F_top(:,2) <= size(indexTopCorr,1),2) =...
            indexTopCorr(F_top(F_top(:,2) <= size(indexTopCorr,1),2),1);
        F_top_shift(F_top(:,3) <= size(indexTopCorr,1),3) =...
            indexTopCorr(F_top(F_top(:,3) <= size(indexTopCorr,1),3),1);
        
%         F_top_shift(F_top(:,1) > size(indexTopCorr,1),1) =...
%             indexTopConc(F_top(F_top(:,1) > size(indexTopCorr,1),1)...
%             - size(indexTopCorr,1),1) + size(V_steelcorr_surf,1);
%         F_top_shift(F_top(:,2) > size(indexTopCorr,1),2) =...
%             indexTopConc(F_top(F_top(:,2) > size(indexTopCorr,1),2)...
%             - size(indexTopCorr,1),1) + size(V_steelcorr_surf,1);
%         F_top_shift(F_top(:,3) > size(indexTopCorr,1),3) =...
%             indexTopConc(F_top(F_top(:,3) > size(indexTopCorr,1),3)...
%             - size(indexTopCorr,1),1) + size(V_steelcorr_surf,1);

        F_top_shift(F_top(:,1) > size(indexTopCorr,1),1) =...
            indexTopConc(F_top(F_top(:,1) > size(indexTopCorr,1),1)...
            - size(indexTopCorr,1),1) + size(Points_corr_surf,1);
        F_top_shift(F_top(:,2) > size(indexTopCorr,1),2) =...
            indexTopConc(F_top(F_top(:,2) > size(indexTopCorr,1),2)...
            - size(indexTopCorr,1),1) + size(Points_corr_surf,1);
        F_top_shift(F_top(:,3) > size(indexTopCorr,1),3) =...
            indexTopConc(F_top(F_top(:,3) > size(indexTopCorr,1),3)...
            - size(indexTopCorr,1),1) + size(Points_corr_surf,1);
        
        F_bot_shift = zeros(size(F_bot,1),3);
        F_bot_shift(F_bot(:,1) <= size(indexBotCorr,1),1) =...
            indexBotCorr(F_bot(F_bot(:,1) <= size(indexBotCorr,1),1),1);
        F_bot_shift(F_bot(:,2) <= size(indexBotCorr,1),2) =...
            indexBotCorr(F_bot(F_bot(:,2) <= size(indexBotCorr,1),2),1);
        F_bot_shift(F_bot(:,3) <= size(indexBotCorr,1),3) =...
            indexBotCorr(F_bot(F_bot(:,3) <= size(indexBotCorr,1),3),1);
        
%         F_bot_shift(F_bot(:,1) > size(indexBotCorr,1),1) =...
%             indexBotConc(F_bot(F_bot(:,1) > size(indexBotCorr,1),1)...
%             - size(indexBotCorr,1),1) + size(V_steelcorr_surf,1);
%         F_bot_shift(F_bot(:,2) > size(indexBotCorr,1),2) =...
%             indexBotConc(F_bot(F_bot(:,2) > size(indexBotCorr,1),2)...
%             - size(indexBotCorr,1),1) + size(V_steelcorr_surf,1);
%         F_bot_shift(F_bot(:,3) > size(indexBotCorr,1),3) =...
%             indexBotConc(F_bot(F_bot(:,3) > size(indexBotCorr,1),3)...
%             - size(indexBotCorr,1),1) + size(V_steelcorr_surf,1);

        F_bot_shift(F_bot(:,1) > size(indexBotCorr,1),1) =...
            indexBotConc(F_bot(F_bot(:,1) > size(indexBotCorr,1),1)...
            - size(indexBotCorr,1),1) + size(Points_corr_surf,1);
        F_bot_shift(F_bot(:,2) > size(indexBotCorr,1),2) =...
            indexBotConc(F_bot(F_bot(:,2) > size(indexBotCorr,1),2)...
            - size(indexBotCorr,1),1) + size(Points_corr_surf,1);
        F_bot_shift(F_bot(:,3) > size(indexBotCorr,1),3) =...
            indexBotConc(F_bot(F_bot(:,3) > size(indexBotCorr,1),3)...
            - size(indexBotCorr,1),1) + size(Points_corr_surf,1);
        
%         F_conc = [F_corr_surf; F_conc_outer + size(V_steelcorr_surf,1); 
%             F_top_shift; F_bot_shift];
        F_conc = [F_corr_surf; F_conc_outer + size(Points_corr_surf,1); 
            F_top_shift; F_bot_shift];