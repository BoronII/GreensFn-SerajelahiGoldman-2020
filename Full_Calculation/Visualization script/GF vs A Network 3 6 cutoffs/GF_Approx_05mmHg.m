% Baran Serajelahi
% Janruary 29th, 2020.
% This script visualizes the tissue oxygen field that is calculatd by
% Greens.cpp

clear all
clear all
close all
close all
%
fs = 12; %font size
fst = 6; %font size for ticks
wid = 4; %line width
NumTicks = 5; %the number of ticks in plot axis

% below I prepare the data to be visualized.
% the data comes as a list
% the first three things listed are the mean, min, max of the oxygen field
% let mxxL, myyL, mzzL be the number of points in the x, y, z directions
% starting from the 4th line of the data file the values of the oxygen
% field are stored
% so there are (mxxL*myyL*mzzL + 3) total values in the data file
% this list moves through the values according to the pattern below
% knowing this will make it possible to understand how the plotting is
% done.
% (1,1,1), (1,1,2)..., (1,1,mzzL),
% (1,2,1), (1,2,2)..., (1,2,mzzL),
% (1,3,1), (1,3,2)..., (1,3,mzzL),
%  ...     ...         ...
% (1,myyL,1), (1,myyL,2)..., (1,myyL,mzzL),
% (2,1,1), (2,1,2)..., (1,1,mzzL),
% (2,2,1), (2,2,2)..., (1,2,mzzL),
% (2,3,1), (2,3,2)..., (1,3,mzzL),
%  ...     ...         ...
% (2,myyL,1), (2,myyL,2)..., (2,myyL,mzzL),
% (3,1,1), (3,1,2)..., (3,1,mzzL),
% (3,2,1), (3,2,2)..., (3,2,mzzL),
% (3,3,1), (3,3,2)..., (3,3,mzzL),
%  ...     ...         ...
%  ...     ...         ...
% (mxxL,1,1), (mxxL,1,2)..., (mxxL,1,mzzL),
% (mxxL,2,1), (mxxL,2,2)..., (mxxL,2,mzzL),
%  ...     ...         ...
%  ...     ...         ...
% (mxxL,myyL,1), (mxxL,myyL,2)..., (mxxL,myyL,mzzL).


% The full GF calculation.
str = sprintf('TissueLevels_NA.dat'); % NA stands for No Approximation cutoff
GF = load(str);
clear str;
GF_mean = GF(1,1);
GF_min 	= GF(2,1);
GF_max 	= GF(3,1);
GF      = GF(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxL = 20;
myyL = 20;
mzzL = 10;
 for i=1:1:mxxL
   for j=1:1:myyL
        GF2(i,j) = GF((6) + mzzL*myyL*(i-1) + (j-1)*mzzL); 
   end
 end


% data for the approximated GF at 198 micron (Low) cuttoff distance
str = sprintf('TissueLevels_AL.dat'); %AL stands for Aproximate Low cutoff
GFAL = load(str); 
clear str;
GFAL_mean = GFAL(1,1);
GFAL_min = GFAL(2,1);
GFAL_max = GFAL(3,1);
GFAL     = GFAL(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxL = 20; %The number of points in the x direction (Can be found/adjusted in simpleNetwork.txt)
myyL = 20; %The number of points in the y direction (Can be found/adjusted in simpleNetwork.txt)
mzzL = 10; %The number of points in the z direction (Can be found/adjusted in simpleNetwork.txt)
 for i=1:1:mxxL
   for j=1:1:myyL
        GFAL2(i,j) = GFAL((6) + mzzL*myyL*(i-1) + (j-1)*mzzL); %This captures the 6th out of the 10 z cross-sections and stores it in a matrix for plotting.
   end
 end
% data for the AL comparison figure
GFLD = abs((GF2 - GFAL2)); 
 
% the approximated GF calculation at 237 micron cuttoff distance. 
str = sprintf('TissueLevels_ALM.dat');% ALM stands for Aproximate Low -Medium cutoff
GFALM = load(str);
clear str;
GFALM_mean = GFALM(1,1);
GFALM_min  = GFALM(2,1);
GFALM_max  = GFALM(3,1);
GFALM      = GFALM(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxM = 20;
myyM = 20;
mzzM = 10;
 for i=1:1:mxxM
   for j=1:1:myyM
        GFALM2(i,j) = GFALM((6) + mzzM*myyM*(i-1) + (j-1)*mzzM); 
   end
 end 
% data for the ALM comparison figure
GFLMD = abs((GF2 - GFALM2)); 

% the approximated GF calculation at 297 micron cuttoff distance. 
str = sprintf('TissueLevels_AM.dat');% AM stands for Aproximate Medium cutoff
GFAM = load(str);
clear str;
GFAM_mean = GFAM(1,1);
GFAM_min  = GFAM(2,1);
GFAM_max  = GFAM(3,1);
GFAM      = GFAM(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxH = 20;
myyH = 20;
mzzH = 10;
 for i=1:1:mxxH
   for j=1:1:myyH
        GFAM2(i,j) = GFAM((6) + mzzH*myyH*(i-1) + (j-1)*mzzH); 
   end
 end
% data for the AM comparison figure
GFMD = abs((GF2 - GFAM2)); 

 % the approximated GF calculation at 396 micron cuttoff distance. 
str = sprintf('TissueLevels_AMH.dat');% AMH stands for Aproximate Medium-High cutoff
GFAMH = load(str);
clear str;
GFAMH_mean 	= GFAMH(1,1);
GFAMH_min 	= GFAMH(2,1);
GFAMH_max 	= GFAMH(3,1);
GFAMH       = GFAMH(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxH = 20;
myyH = 20;
mzzH = 10;
 for i=1:1:mxxH
   for j=1:1:myyH
        GFAMH2(i,j) = GFAMH((6) + mzzH*myyH*(i-1) + (j-1)*mzzH);    
   end
 end
% data for the AMH comparison figure
GFMHD =  abs((GF2 - GFAMH2)); 
 
% the approximated GF calculation at 495 micron cuttoff distance. 
str = sprintf('TissueLevels_AH.dat');% AH stands for Aproximate High cutoff
GFAH = load(str);
clear str;
GFAH_mean = GFAH(1,1);
GFAH_min  = GFAH(2,1);
GFAH_max  = GFAH(3,1);
GFAH      = GFAH(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxH = 20;
myyH = 20;
mzzH = 10;
 for i=1:1:mxxH
   for j=1:1:myyH
        GFAH2(i,j) = GFAH((6) + mzzH*myyH*(i-1) + (j-1)*mzzH); 
   end
 end
% data for the AH comparison figure
GFHD =  abs((GF2 - GFAH2)); 

% the approximated GF calculation at 59418 micron cuttoff distance. 
str = sprintf('TissueLevels_AVH.dat');% AVH stands for Aproximate Very High cutoff
GFAVH = load(str);
clear str;
GFAVH_mean = GFAVH(1,1);
GFAVH_min  = GFAVH(2,1);
GFAVH_max  = GFAVH(3,1);
GFAVH      = GFAVH(4:end,1); %Drop the first 3 entries (mean, min and, max) and just keep the tissue data
mxxH = 20;
myyH = 20;
mzzH = 10;
 for i=1:1:mxxH
   for j=1:1:myyH
        GFAVH2(i,j) = GFAVH((6) + mzzH*myyH*(i-1) + (j-1)*mzzH); 
   end
 end
% data for the VHA comparison figure
GFVHD =  abs((GF2 - GFAVH2)); 

% Now lets plot the results in a 3*3 panel figure.
figure('rend','painters','pos',[1 1 500 500]);
hold on
% Panel A
subplot (6,3,1)
contourf(GFAL2,'LineStyle','none')
box off
title('Cutoff 198 microns, 1.15 s/i, 38 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel B
subplot (6,3,2)
contourf(GF2,'LineStyle','none')
box off
title('Full, 1.07 s/i, 27 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -6;
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
% txt1 			= '0'; %tick labels for x axis
% text(1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(2*3.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(2*9.25, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(2*13.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(2*18 - 1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0'; % tick labels for y axis
% text(-2, 1, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(-5.5, 2*5, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(-4, 2*10.25, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(-5.5, 2*15.25, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(-4, 2*19 - 1, txt1, 'Fontsize', fs);

% Panel C
subplot (6,3,3)
contourf(GFLD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel D
subplot (6,3,4)
contourf(GFALM2,'LineStyle','none')
box off
title('Cutoff 237 microns, 1.18 s/i, 59 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;


% Panel E
subplot (6,3,5)
contourf(GF2,'LineStyle','none')
box off
title('Full')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel F
subplot (6,3,6)
contourf(GFLMD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel G
subplot (6,3,7)
contourf(GFAM2,'LineStyle','none')
box off
title('Cutoff 297 microns, 1.08 s/i, 49 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -6;
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
% txt1 			= '0'; %tick labels for x axis
% text(1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(2*3.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(2*9.25, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(2*13.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(2*18 - 1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0'; % tick labels for y axis
% text(-2, 1, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(-5.5, 2*5, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(-4, 2*10.25, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(-5.5, 2*15.25, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(-4, 2*19 - 1, txt1, 'Fontsize', fs);

% Panel H
subplot (6,3,8)
contourf(GF2,'LineStyle','none')
box off
title('Full')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel I
subplot (6,3,9)
contourf(GFMD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max Residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel K
subplot (6,3,10)
contourf(GFAMH2,'LineStyle','none')
box off
title('Cutoff 396 microns, 1.05 s/i, 32 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel L
subplot (6,3,11)
contourf(GF2,'LineStyle','none')
box off
title('Full')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -6;
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
% txt1 			= '0'; %tick labels for x axis
% text(1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(2*3.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(2*9.25, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(2*13.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(2*18 - 1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0'; % tick labels for y axis
% text(-2, 1, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(-5.5, 2*5, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(-4, 2*10.25, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(-5.5, 2*15.25, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(-4, 2*19 - 1, txt1, 'Fontsize', fs);

% Panel M
subplot (6,3,12)
contourf(GFMHD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel N
subplot (6,3,13)
contourf(GFAH2,'LineStyle','none')
box off
title('Cutoff 495 microns, 1.08 s/i, 26 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel O
subplot (6,3,14)
contourf(GF2,'LineStyle','none')
box off
title('Full')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -6;
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
% txt1 			= '0'; %tick labels for x axis
% text(1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(2*3.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(2*9.25, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(2*13.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(2*18 - 1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0'; % tick labels for y axis
% text(-2, 1, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(-5.5, 2*5, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(-4, 2*10.25, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(-5.5, 2*15.25, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(-4, 2*19 - 1, txt1, 'Fontsize', fs);

% Panel P
subplot (6,3,15)
contourf(GFHD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel Q
subplot (6,3,16)
contourf(GFAVH2,'LineStyle','none')
box off
title('Cutoff 59418 microns, 1.05 s/i, 27 i')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;

% Panel R
subplot (6,3,17)
contourf(GF2,'LineStyle','none')
box off
title('Full')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Pressure (mmHg)';
c.Box = 'off';
caxis([0 50]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -6;
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
% txt1 			= '0'; %tick labels for x axis
% text(1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(2*3.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(2*9.25, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(2*13.75, -3, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(2*18 - 1, -3, txt1, 'Fontsize', fs);
% txt1 			= '0'; % tick labels for y axis
% text(-2, 1, txt1, 'Fontsize', fs);
% txt1 			= '0.005';
% text(-5.5, 2*5, txt1, 'Fontsize', fs);
% txt1 			= '0.01';
% text(-4, 2*10.25, txt1, 'Fontsize', fs);
% txt1 			= '0.015';
% text(-5.5, 2*15.25, txt1, 'Fontsize', fs);
% txt1 			= '0.02';
% text(-4, 2*19 - 1, txt1, 'Fontsize', fs);

% Panel S
subplot (6,3,18)
contourf(GFVHD,'LineStyle','none')
box off
title('Difference')
% xl = xlabel('x (cm)');
% ylabel('y (cm)')
set(gca, 'FontSize', fs);
colormap(jet)
colorbar
c = colorbar;
c.Label.String = 'Max residue (mmHg)';
c.Box = 'off';
caxis([0 20]);
set(gca, 'FontSize', fs);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = wid;
xl.Position(2) = -9;




str = sprintf('GF_Approx.tif');
saveas(gcf, str);
clear str;







