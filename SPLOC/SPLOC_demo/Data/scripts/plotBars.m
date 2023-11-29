%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    plotBars.m                           %%%                                   %%%
%%%  Purpose: Generate bar plots from SPLOC runs   %%%    BioMolecular Physics Group     %%%
%%%           under different conditions           %%%   University of North Carolina    %%%
%%%  Created: 07-23-2023                           %%%           at Charlotte            %%%
%%%  Author:  Tyler J Grear                        %%%                                   %%%
%%----------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp1 = [0 15 0];
grp2 = [0 15 0];
grp3 = [2 2 11];
grp4 = [3 2 10];
X = vertcat(grp1,grp2,grp3,grp4);
high = max(max(X));
Y = categorical({'Case 1','Case 2','Case 3','Case 4'});
Y = reordercats(Y,{'Case 1','Case 2','Case 3','Case 4'});
%==========================================================================================%

figure(5040)
b = bar(Y,X,'stacked');
hold on
%title("A Titely Title is Titlrific");
%legend("d-modes","u-modes","i-modes")
ylabel("Number of modes");
ylim([0 (high + 1.2)]);

b(1).FaceColor = 'flat';
b(1).CData(1,:) = [1 0 0];
b(2).FaceColor = 'flat';
b(2).CData(1,:) = [1 1 0];
b(3).FaceColor = 'flat';
b(3).CData(1,:) = [0 0 1];

b(1).FaceColor = 'flat';
b(1).CData(2,:) = [1 0 0];
b(2).FaceColor = 'flat';
b(2).CData(2,:) = [1 1 0];
b(3).FaceColor = 'flat';
b(3).CData(2,:) = [0 0 1];

b(1).FaceColor = 'flat';
b(1).CData(3,:) = [1 0 0];
b(2).FaceColor = 'flat';
b(2).CData(3,:) = [1 1 0];
b(3).FaceColor = 'flat';
b(3).CData(3,:) = [0 0 1];

b(1).FaceColor = 'flat';
b(1).CData(4,:) = [1 0 0];
b(2).FaceColor = 'flat';
b(2).CData(4,:) = [1 1 0];
b(3).FaceColor = 'flat';
b(3).CData(4,:) = [0 0 1];

legend("d-modes","u-modes","i-modes",'location','southoutside','orientation','horizontal')

%{
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
%}

hold off

set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]); delete("temp"+mode+".png");
exportgraphics(gcf,"4-barPlot_stacked_2.png","Resolution",300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%