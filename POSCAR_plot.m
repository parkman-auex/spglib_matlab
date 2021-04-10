
%% 
function fig  = POSCAR_plot(Rm,sites,Atom_name,Atom_num)
% 
Standradio =30;
%
if nargin <1
    POSCAR_read;
end
fig = figure();
hold on, axis equal
axis image on
grid on;
% supercell first
    
% 
plotBox(Rm);

% 
 elements_table=element_table();
 %[~,sites_num]=size(sites);
 
 tempseq=0;
 for i=1:length(Atom_name)
     for j=1:Atom_num(i)
        %plotAtoms(position2, [244, 13, 100]/255, 40);
        tempseq    = tempseq+1;
        position   = [sites(tempseq).rc1,sites(tempseq).rc2,sites(tempseq).rc3];
        % disp(position );
        sites_name = Atom_name(i);
        temprows = ( elements_table.element_name==sites_name);
        %disp(temprows);
        markercolor = table2array(elements_table(temprows,{'rgb1','rgb2','rgb3'}));
        markersize =  table2array(elements_table(temprows,{'element_radius'}))*Standradio;
        plotAtoms(Rm,position, markercolor, markersize); 
     end
 end
 
titlename = "";
for i =1:length(Atom_name)
    titlename=titlename+Atom_name(i)+Atom_num(i);
end

title(char(titlename));
view(60,10);
axis equal;
end

%% function
function [a,b,c,alpha,beta,gamma]=Rm2abc(Rm)
    a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    a=norm(a1);b=norm(a2);c=norm(a3);
end
function plotBox(Rm)
    %fig = figure();
    import vasplib_tool_outer.*;
    %box on;
    a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    % 8site
    vertex = [0, 0, 0;...
              [a1];...
              [a2];...
              [a1+a2];...
              [a3];...
              [a1+a3];...
              [a2+a3];...
              [a1+a2+a3]];
    % 12
    ax = gca;
    arrow3(vertex(1,:),vertex(2,:)*0.2);
    arrow3(vertex(1,:),vertex(3,:)*0.2);
    arrow3(vertex(1,:),vertex(5,:)*0.2);
    plotLine(ax,vertex(1,:), vertex(2,:));
    plotLine(ax,vertex(1,:), vertex(3,:));
    plotLine(ax,vertex(2,:), vertex(4,:));
    plotLine(ax,vertex(3,:), vertex(4,:));
    plotLine(ax,vertex(5,:), vertex(6,:));
    plotLine(ax,vertex(5,:), vertex(7,:));
    plotLine(ax,vertex(6,:), vertex(8,:));
    plotLine(ax,vertex(7,:), vertex(8,:));
    plotLine(ax,vertex(1,:), vertex(5,:));
    plotLine(ax,vertex(2,:), vertex(6,:));
    plotLine(ax,vertex(3,:), vertex(7,:));
    plotLine(ax,vertex(4,:), vertex(8,:));
    % 
%     Rm_=-Rm;
%     [x_min, x_max, y_min, y_max, z_min, z_max] = deal(min(Rm_(:,1)), max(Rm(:,1)), ...
%                                                       min(Rm_(:,2)), max(Rm(:,2)), ...
%                                                       min(Rm_(:,3)), max(Rm(:,3)));
%     x_len = (x_max - x_min)/2;
%     y_len = (y_max - y_min)/2;
%     z_len = (z_max - z_min)/2;
%    
%     axis([x_min-x_len x_max+x_len y_min-y_len y_max+y_len z_min-z_len z_max+z_len]);
     
end

function plotAtoms(Rm,position, markercolor, markersize)
    % a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    real_position=position*Rm;
    cur_point = real_position;
    ax =gca;
    plot3(ax,cur_point(:,1), cur_point(:,2), cur_point(:,3), 'ok', 'linewidth', 1.5, 'markersize', markersize, 'markerfacecolor', markercolor)

end

function plotLine(ax,x1, x2)
%   
  
    plot3(ax,[x1(1) x2(1)], [x1(2), x2(2)], [x1(3), x2(3)], 'k', 'linewidth', 1.3)
end
%%
function elements_table=element_table()
%% Initialize variables.
filename = '/Users/parkman/Documents/MATLAB/sourcefile/elements.dat';
delimiter = ' ';

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
elements_table = table(dataArray{1:end-1}, 'VariableNames', {'element_name','element_value','element_radius','var1','var2','rgb1','rgb2','rgb3'});


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
%%
end