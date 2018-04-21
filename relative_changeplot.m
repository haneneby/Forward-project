%%%%%% relative change file plot results %%%%%%5
close all
clear all
M = csvread ('/Users/HKhanene/Documents/MATLAB/back_up/aRC/RC_newhigh.csv');
M = csvread ('/Users/HKhanene/Documents/MATLAB/back_up/aRC/RC_newlow.csv');

A=size(M,2)
Baseline=M(:,1:1);% lesion move by D
for i= 2:A
 B=  M(:,i:i);
 C= real(B)-real(Baseline);
 C= C./real(Baseline);
 size(C)
 if (i==2)
     max_results= max(C);
     results= C;
 else
    max_results= cat(2,max_results,max(C));
    results= cat(2,results,C);

 end
 
end
'results'
size(results)
 D = [0.5 1 1.5 2 2.5 3 3.5 4.5 5 6 ]
%  D = [3 4 5 6 8 10 11];% lesion move by D
% D = [15 20 25 30 35 40 45 50];% lesion translate by D
 %rr= cat(2,max_results(:,8),max_results(:,1:7));
 figure;
%  plot( D,real (rr))% lession size increase by D
 plot( D,real (max_results));% lesion move by D
% xlabel('Lesion  size increases by a factor D '); % lession size increase by D
  xlabel('Lesion  is at a distance D to source 2'); % lession size increase by D

 %xlabel('Lesion has a depth D in mm ');% lesion move by D
 ylabel('Relative change');
marker=['-+','-o','-*','..','-x','s','d','-^','v','>','<','p','h'];
 figure;
 hold all
 for k = 1:length(D) 
 plot(real (results(:,k)),marker(k));
 end
legend('x0.5','x1','x1.5','x2', 'x2.5','x3','x3.5','x4.5','x5','x6');

%  D = [9 11 16 21 26 31 36 41];% lesion translate by D
% 
%  xlabel('Detector index: S1 only then S2 only ');
%  ylabel('Relative change');
% %  legend( '3mm', '4mm','5mm','6mm','8mm','10mm','11mm'); % lesion move by D
% % legend('15mm','20mm','25mm','30mm','35mm','40mm','45mm','50mm '); % lesion translate by D


%Baseline=M(:,7:7);% lession size increase by D

% R1= M(:,1:6); % lession size increase by D
% R2= M(:,8:11);% lession size increase by D
% R21=cat(2,R2,R1);% lession size increase by D
