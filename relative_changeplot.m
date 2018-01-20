%%%%%% relative change file plot results %%%%%%5
close all
clear all
M = csvread ('/Users/HKhanene/Documents/MATLAB/parser_test/bcRC/RC.csv');
A=size(M,2)
Baseline=M(:,1:1);% lesion move by D
% Baseline=M(:,7:7);% lession size increase by D

% R1= M(:,1:6); % lession size increase by D
% R2= M(:,8:10);% lession size increase by D
% R21=cat(2,R2,R1);% lession size increase by D
% j=1;
for i= 2:A%1:A-1
 B=  M(:,i:i);%R21(:,i:i);
 C= real(Baseline)-real(B);
 C= C./real(Baseline);
 size(C)
 if (i==2) %(i==1)
     max_results= max(C);
     results= C;

     j=2;
 else
    max_results= cat(2,max_results,max(C));
    results= cat(2,results,C);

 end
 
end
'results'
size(results)
% D = [0.5 1 1.5 2 2.5 3 3.5 4.5  6 ]
%  D = [3 4 5 6 8 10 11];% lesion move by D
 D = [15 20 25 30 35 40 45 50];% lesion translate by D
 rr= cat(2,max_results(:,8),max_results(:,1:7));
 figure;
 plot( D,real (rr))% lession size increase by D
 %plot( D,real (max_results));% lesion move by D
% xlabel('Lesion  size increases by a factor D '); % lession size increase by D
  xlabel('Lesion  is at a distance D to source 2'); % lession size increase by D

 %xlabel('Lesion has a depth D in mm ');% lesion move by D
 ylabel('Relative change');
 
 figure;
 plot(real (results(:,8)),'-x');
 hold on
 plot(real (results(:,1)),'- +');
 hold on
 plot(real (results(:,2)),' *');
  hold on
 plot(real (results(:,3)),'-.');
  hold on
 plot(real (results(:,4)),' - ');
  hold on
 plot(real (results(:,5)),' -.*');
  hold on
 plot(real (results(:,6)),'--');
  hold on
 plot(real (results(:,7)),'m .');
%  hold on
%  plot(real (results(:,8)),'m.');
%  hold on
%  plot(real (results(:,9)),'y p-');
%  hold on
% plot(real (results(:,10)),'g*');
 D = [9 11 16 21 26 31 36 41];% lesion translate by D

 xlabel('Detector index: S1 only then S2 only ');
 ylabel('Relative change');
%  legend('x0.5','x1','x1.5','x2', 'x2.5','x3','x3.5','x4.5','x6');
%  legend( '3mm', '4mm','5mm','6mm','8mm','10mm','11mm'); % lesion move by D
 legend('15mm','20mm','25mm','30mm','35mm','40mm','45mm','50mm '); % lesion translate by D

