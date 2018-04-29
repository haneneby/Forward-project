%-----------------------------------------------------------------------%
%----- data to .geo file --------------------------------------------%
%-----------------------------------------------------------------------%

clc
close all
clear

%-----------------------------------------------------------------------%
% define border
% lc   = 3;                        % mesh size
% Xc   = [ -25,25,-25,25];         % X - coordinates of corners
% Zc   = [ 0,0,-50,-50];           % Z - coordinates of corners

rad =60;

for k=2351:2400

Xc   = [ 0,60,0, 60]; %15, 17+45.* rand(1)       % X - coordinates of corners
Zc   = [ 0,0,60, 60];           % ,10+50.* rand(1) y - coordinates of corners

A(1:4,1:5) = 0;
for i =1:4
   A(i,1)=i ; 
   A(i,2)=Xc(i);   
   A(i,3)=Zc(i);                      % Y-coordinate (X-Z plane)
   A(i,4)=0;
   A(i,5)=1;
end

%----------------------------------------------------------------------%
%- lines that join the border points

Prev = [ 1;1;3;4];   % previous corner node
Next = [ 2;3;4;2];   % next corner node

L(1:4,1:3) = 0;
for i = 1:4
    L(i,1) = i;
    L(i,2) = Prev(i);
    L(i,3) = Next(i);
end
  
%-----------------------------------------------------------------------%
% internal points where finer mesh is needed

lc1   = 0.3;                     % mesh size
X    = [ -2,2,1,-1];             % X - coordinates of internal box
Z    = [ -1,-1,-3,-3];           % Z - coordinates of internal box

B(1:4,1:5) = 0;
for i =1:4
   B(i,1)=i+4 ; 
   B(i,2)=X(i);   
   B(i,3)=0;                      % Y-coordinate (X-Z plane)
   B(i,4)=Z(i);
   B(i,5)=lc1;
end
 % internal blod where finer surface is needed
   
    for i =1:1
     a   = 2;%2+2*rand(1);  % ellipses large diameter [2 10] was  2+8
     b   = 2;%2+(a-2)*rand(1);   % ellipses small diameter [0 a]
    if (i==1)
     Y= 25+ a + (0.2*rad-a-b).*rand(1);   %0.4*rad y of the center of ellipse 
    else
      Y= a + (0.*rad-a-b).*rand(1)+ D(i-1,3)./2;%0.4*rad  if upper shape is thinner
    end
     %if (Y > 0.2*rad)
%         if (i==1)
%             X= 0.46*rad-a +(0.2.*rad+a-10).*rand(1); % x of the center of ellipse 
%         else
%             X= 0.46*rad-a +(0.2.*rad+a-10).*rand(1)-12;
%         end
%         
    %else
        if (i==1)
            X= a +((rad-2*a)).*rand(1);
        
        else
          X= 0.46*rad-a +(0.2.*rad+a-10).*rand(1)+12;
        end

    %end
       D(i,1)=i+4; 
       D(i,2)=fix (X)+(i-1)*a;   
       D(i,3)=fix(Y);                      % Y-coordinate (X-Z plane)
       D(i,4)=0;
       D(i,5)=a;
       D(i,6)=b;
       D(i,7)=0;
       %D(i,8)='2*Pi';
    end
% GMsh geometry file is stored as .txt file
filename = '/Users/HKhanene/Documents/ML_dataset/phanshapes/%d';
    filename = sprintf(filename,k);
    filename= strcat(filename,'.geo');
    fileID = fopen(filename,'w');
for i=1:4
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',A(i,:));
end
for i = 1:4
fprintf(fileID,'Line(%d) = {%d,%d};\n',L(i,:));
end
fprintf(fileID,'Line Loop(2) = {1,-4,-3,-2};\n');
fprintf(fileID,'Plane Surface(1) = {2};\n');
fprintf(fileID,'SetFactory("OpenCASCADE");\n');
fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(1,:));
%fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(2,:));
%fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(3,:));
fprintf(fileID,'Line Loop(3) = {5};\n');
fprintf(fileID,'Plane Surface(2) = {3};\n');
%fprintf(fileID,'Line Loop(4) = {6};\n');
%fprintf(fileID,'Plane Surface(3) = {4};\n');
% fprintf(fileID,'Line Loop(5) = {7};\n');
% fprintf(fileID,'Plane Surface(4) = {5};\n');
fclose(fileID);
% type geofile.txt
% uiopen('C:\ROOT GEOMETRY\RootSys Surface\geofile.txt',1)
end
%-----------------------------------------------------------------