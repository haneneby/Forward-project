%%%%%
% verify mesh read save figure and result in a files and csv

%%%

clear all
close all

% myDir = uigetdir; %gets directory
% myFiles = dir(fullfile(myDir,'*.txt'); %gets all txt files in struct
% for k = 1:length(myFiles)
%   baseFileName = myFiles(k).name;
%   fullFileName = fullfile(myDir, baseFileName);
%   fprintf(1, 'Now reading %s\n', fullFileName);
%   num = importdata(fullFileName);   %or readtable
%   % all of your actions for filtering and plotting go here
% end

%file = [breast_dim.msh triangle_2.msh]
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/gmsh_test/train/triangle_2.msh','gmsh');
mesh.Display

results= zeros([128 1]);
ll=0;
% source detectore 
Q(1,:)=[0 0];%[6 -27];
Q(2,:)=[60 0];
Q(4);
d= sqrt((Q(1)-Q(2))^2 +(Q(3)-Q(4))^2);
d1=[16.33 45];
t1=d1(1)/d;%Let the ratio of distances, 
t2=d1(2)/d;%Let the ratio of distances, 
M(1,:)= [(1-t1)*Q(1)+t1*Q(2) (1-t1)*Q(3)+t1*Q(4)];
M(128,:)= [(1-t2)*Q(1)+t2*Q(2) (1-t2)*Q(3)+t2*Q(4)];
V(1,:)=M(1,:);
V(2,:)=M(128,:);
dis=28.67/128;
di= t1;
for i=2:128
ti=di/28.67;%Let the ratio of distances, 

M(i,:)= [(1-ti)*V(1)+ti*V(2) (1-ti)*V(3)+ti*V(4)];
di= di+dis;

end

% system matrix
n = mesh.NodeCount;
a= 0.002;
b= 0.012;
a1= 0.002;
b1= 0.012;
firstloop=0;
for freq=[352.94 375 400 434.78]
%create different meshes with various absorption and scattering
%coefficients
for j=0:100
mua = (b-a).* rand(n,1)*0.01 +a; % absorption coefficient
mus = (b1-a1).* rand(n,1)*1+a1; % scattering coefficient
ref = rand(n,1)*1.4; % refractive index
%
% Map to mesh basis
%figure;
figure('visible','off');
mesh.Display(mua);
set(gca,'XTick',[]) % Remove the ticks in the x axis!
set(gca,'YTick',[]) % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
c = colorbar;
c.Visible = 'off';
filename = '/Users/HKhanene/Documents/MATLAB/train_set/y_train/meshabsorption%d-%d';
filename = sprintf(filename,fix(freq),401+j+ll);
saveas(gcf,filename,'png')
%figure; mesh.Display(mus);
Y= zeros([128 1]);
for i=1:2
mesh.SetQM(Q(i,:),M);
%display source /detectors
 %hold on
%  plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
%  plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');
%  %bb=mesh.BoundingBox()

nnd = mesh.NodeCount;
ref_bkg = 1.4;
ref = ones(n,1) * ref_bkg;

% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2, ref);

% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,freq);
Phi = K\qvec;
Y(:,i) = mvec.' * Phi;
end
size(abs(Y));

% flatten the Y vector for ML processing
% MESVEC= cat(1,Y(:,1),Y(:,2));
% superposition of two light path lengths of illuminated light from eLED #1 and eLED #2.
MESVEC= sum(Y,2);

size(abs(MESVEC));

logY = log(Y);
lnamp = real(logY);
phase = imag(logY);
% % Display sinogram
% figure
% subplot(1,2,1); 
% imagesc(lnamp);
% title('logarithmic amplitude');
% xlabel('source index q');
% ylabel('detector index m');
% axis equal tight;
% colorbar
% subplot(1,2,2); 
% imagesc(phase);
% title('phase');
% xlabel('source index q');
% ylabel('detector index m');
% axis equal tight;
% colorbar

% size(abs(Y))
% % Write solver results to file
% data = reshape(Y,[128,2]);
% file = sprintf('/Users/HKhanene/Documents/MATLAB/train_set/x_train/measurement%d', j),
% file= strcat(file,'.dat')
% toastWriteVector(file,Y);
% if(j==0 && firstloop==0)
%     csvwrite('/Users/HKhanene/Documents/MATLAB/train_set/breast.csv',full(MESVEC),0,0);
%     j
% else
%     datacsv = csvread('breast.csv');
%     newdata = [datacsv full(MESVEC)];
%     csvwrite('/Users/HKhanene/Documents/MATLAB/train_set/breast.csv',newdata);
%     j
% end
 %dlmwrite('breast.csv',full(MESVEC),'-append');
 if(j==0 && firstloop==0)
    results= MESVEC;
 else
    results= cat(2,results,MESVEC);
 end

end
firstloop=1;
size(results);
ll=ll+100

end
csvwrite('/Users/HKhanene/Documents/MATLAB/train_set/breast1.csv',full(results),0,0);

datacsv = csvread('/Users/HKhanene/Documents/MATLAB/train_set/breast.csv');
newdata = [datacsv full(MESVEC)];
csvwrite('/Users/HKhanene/Documents/MATLAB/train_set/breast2.csv',newdata);


