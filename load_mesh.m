clear all
close all
a= 0.002;
b= 0.012;
a1= 0.002;
b1= 0.012;
count= 0;
mm=0;
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.msh')); %gets all txt files in struct
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  num = importdata(fullFileName);   %or readtable
  % all of your actions for filtering and plotting go here
  
  mesh = toastMesh(fullFileName,'gmsh');
mesh.Display

% system matrix
n = mesh.NodeCount;
mua = (b-a).* rand(n,1)*0.01 +a; % absorption coefficient
mus = (b1-a1).* rand(n,1)*1+a1; % scattering coefficient
ref = rand(n,1)*1.4; % refractive index
%
% Map to mesh basis
figure;
mesh.Display(mua);
set(gca,'XTick',[]) % Remove the ticks in the x axis!
set(gca,'YTick',[]) % Remove the ticks in the y axis
set(gca,'Position',[0 0 1.28 1.28]) % Make the axes occupy the hole figure
c = colorbar;
c.Visible = 'off';

filename = '/Users/HKhanene/Documents/MATLAB/train_set/y_test/meshabs%d-%d';
filename = sprintf(filename,k,k);
saveas(gcf,filename,'png')

results= zeros([128 1]);
ll=0;
% source detectore 
% if mm<3 
%     Q(1,:)=[6 -27];
% else
    Q(1,:)=[0 0];%[6 -27];
% end
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


mesh.SetQM(Q,M);
%display source /detectors
 hold on
 plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
 plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');
 %bb=mesh.BoundingBox()

mm= mm+1
end
