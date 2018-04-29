
%%%%%%
% create different .geo file with diff shape
%%%%%

%-----------------------------------------------------------------------%
%----- data to .geo file --------------------------------------------%
%-----------------------------------------------------------------------%

clc
close all
clear

%-----------------------------------------------------------------------%
rad =60;
for k=622:650
   

    %-----------------------------------------------------------------------%
    % internal blod where finer surface is needed
   
    for i =1:3
     a   = 10*rand(1);  % ellipses large diameter [0 10]
     b   = a*rand(1);   % ellipses small diameter [0 a]

     Y= a + (0.4*rad-a-b).*rand(1);   % y of the center of ellipse 
    if (Y > 0.2*rad)
            X= 0.46*rad-a +(0.2.*rad+a).*rand(1); % x of the center of ellipse 
    else
        X= 10+a +(rad-2*a).*rand(1);
    end
       D(i,1)=i+5; 
       D(i,2)=fix (X);   
       D(i,3)=fix(Y);                      % Y-coordinate (X-Z plane)
       D(i,4)=0;
       D(i,5)=a;
       D(i,6)=b;
       D(i,7)=0;
       %D(i,8)='2*Pi';
    end
    % GMsh geometry file is stored as .txt file
    filename = '/Users/HKhanene/Documents/MATLAB/dataset/%d';
    filename = sprintf(filename,k);
    filename= strcat(filename,'.geo');
    fileID = fopen(filename,'w');
    %fprintf(fileID, 'rad = DefineNumber[ 60, Name "Parameters/rad" ];\n);
    fprintf(fileID,'Point(3)={%d,0,0,1};\n',rad);
    fprintf(fileID,'Point(4)={%d,%d,0,1};\n',0.5*rad,0.6*rad);
    fprintf(fileID,'Point(5)={0,0,0,1};\n');
    fprintf(fileID,'Point(6)={0,0,0,1};\n');
    fprintf(fileID,'BSpline(1) = {6, 4, 3};\n');
    fprintf(fileID,'Line(2) = {6, 3};\n');
    fprintf(fileID,'Line Loop(1) = {1,-2};\n');
    fprintf(fileID,'Plane Surface(1) = {1};\n');
    fprintf(fileID,'SetFactory("OpenCASCADE");\n');
    fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(1,:));
    fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(2,:));
    fprintf(fileID,'Ellipse(%d) = {%d,%d,%d,%d,%d,%d,2*Pi};\n',D(3,:));

    fprintf(fileID,'Line Loop(3) = {6};\n');
    fprintf(fileID,'Plane Surface(2) = {3};\n');
    fprintf(fileID,'Line Loop(4) = {7};\n');
    fprintf(fileID,'Plane Surface(3) = {4};\n');
    fprintf(fileID,'Line Loop(5) = {8};\n');
    fprintf(fileID,'Plane Surface(4) = {5};\n');
    fclose(fileID);
    %type filename;
end
%-----------------------------------------------------------------