
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
% define border
X1=100.* rand(100,1)
X2=100.* rand(100,1)
Y1=60.* rand(100,1);
Y2=60.* rand(100,1);

for k=1:10
    lc   = 3;                        % mesh size
    Xc   = [ 0,60,X1(k),X2(k)];        % X - coordinates of corners
    Zc   = [ 0,0,Y1(k),Y2(k)];           % y - coordinates of corners

    A(1:4,1:5) = 0;
    for i =1:4
       A(i,1)=i ; 
       A(i,2)=Xc(i);   
       A(i,3)=Zc(i);                      % Y-coordinate (X-Z plane)
       A(i,4)=0;
       A(i,5)=lc;
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
    X    = [ 28,32,30,29];             % X - coordinates of internal box
    Z    = [ 1,1,3,3];           % Y - coordinates of internal box

    B(1:4,1:5) = 0;
    for i =1:4
       B(i,1)=i+4 ; 
       B(i,2)=X(i);   
       B(i,3)=Z(i);                      % Y-coordinate (X-Z plane)
       B(i,4)=0;
       B(i,5)=lc1;
    end
    %-----------------------------------------------------------------------%
    % internal blod where finer surface is needed
    rad1   = [ 2,5,4,6];                     % mesh size
    rad2   = [ 1,3,2,5];
    X    =10+ abs((((-X1(k)+X2(k)) )/ (rad1(k) *rad2(k))).* rand(4,1))%[ 28,32,13,25]; %e(X1(k)-X2(k)).* rand(4,1)- Y2(k);             % X - coordinates of internal box
    Z    =10+abs((((-Y1(k)+Y2(k)) )/ (rad1(k)* rad2(k))).* rand(4,1))%[ -6,-24,-8,-23]; %(Y1(k)+Y2(k)).* rand(4,1)+ X1(k)% (Y1(k)-Y2(k)).* rand(4,1)- Y2(k);           % Y - coordinates of internal box
     
    D(1:4,1:6) = 0;
    for i =1:4
       D(i,1)=i+4 ; 
       D(i,2)=X(i);   
       D(i,3)=Z(i);                      % Y-coordinate (X-Z plane)
       D(i,4)=0;
       D(i,5)=rad1(i);
       D(i,6)=rad2(i);
    end
    % GMsh geometry file is stored as .txt file
    filename = '/Users/HKhanene/Documents/MATLAB/dataset/geo%d';
    filename = sprintf(filename,k);
    filename= strcat(filename,'.geo');
    fileID = fopen(filename,'w');
    for i=1:4
    fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',A(i,:));
    end
    for i = 1:4
    fprintf(fileID,'Line(1) = {1,2};\n');
    end
    fprintf(fileID,'Line Loop(5) = {1,-4,-3,-2};\n');
    %fprintf(fileID,'Line Loop(5) = {1,-2};\n');
    fprintf(fileID,'Plane Surface(1) = {5};\n');
    for i=1:4
    fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',B(i,:));
    end
    for i=1:4
    fprintf(fileID,'Point{%d} In Surface{6};\n',B(i,1));
    end
    fprintf(fileID,'SetFactory("OpenCASCADE");\n');
    for i=2:3
    fprintf(fileID,'Disk(%d)={%d,%d,%d,%d,%d};\n',D(i,:));
    end
    fprintf(fileID,'Line Loop(8) = {4};\n');
    fprintf(fileID,'Surface(8) = {8};\n');
%     fprintf(fileID,'//+\n');
     fprintf(fileID,'Line Loop(10) = {3};\n');
%     fprintf(fileID,'//+\n');
    fprintf(fileID,'Surface(9) = {10};\n');
%    

%     %for i=1:4
%     i=7;
%     j=7+i;
%     fprintf(fileID,'Line Loop(10) = {9};\n');

    %fprintf(fileID,'Line Loop(%i) = {%j};\n');
    %en
%     for i=1:4
%     j=7+i;
%     fprintf(fileID,'Plane Surface(%j) = {%d};\n',D(i,1));
%     end
  
    fclose(fileID);
    %type filename;
    %uiopen('C:\ROOT GEOMETRY\RootSys Surface\geofile2.txt',1)
end
%-----------------------------------------------------------------