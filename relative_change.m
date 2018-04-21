%%%%%% relative change file%%%%%%5
close all
clear all

% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.002; % background absorption [1/mm]
mus_bkg = 2.5;    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion
w= .00069;        % wave length
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.05;      % scattering coef of medium index range uapper limit
tol = 1e-6;     % tolerance of bicgstab where the default is  1e-6.
maxit= 50 ;     % max number of iteration, default is min (n,20) where n is the size of vector b in bicgstab(A,b)
f = cm/w
loop=0;
results= zeros([128 1]);
bx=128; by=128;

% load mesh directory
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.msh')); %gets all txt files in struct
myFiles = natsortfiles({myFiles.name});
%myFiles= orderfields(myFiles);
for k = 1:length(myFiles)
  baseFileName = myFiles(k)%.name;
  fullFileName = fullfile(myDir,'/', baseFileName)
  %fprintf(1, 'Now reading %s\n', char(fullFileName));
  %num = importdata(fullFileName);  %or readtable
  % all of your actions for filtering and plotting go here
 % load mesh% load mesh

  mesh = toastMesh(char(fullFileName),'gmsh');
  ne = mesh.ElementCount;
  nv = mesh.NodeCount;
  regidx = mesh.Region;
  regno = unique(regidx);
  [vtx,idx,eltp] = mesh.Data;

% define source 
 Q(1,:)=[0 0];%[6 -27];
 Q(2,:)=[60 0];
%define detector locations
d= sqrt((Q(1)-Q(2))^2 +(Q(3)-Q(4))^2);
d1=[16.33 45];
t1=d1(1)/d;% the ratio of distances, 
t2=d1(2)/d;% the ratio of distances, 
M(1,:)= [(1-t1)*Q(1)+t1*Q(2) (1-t1)*Q(3)+t1*Q(4)];
M(128,:)= [(1-t2)*Q(1)+t2*Q(2) (1-t2)*Q(3)+t2*Q(4)];
V(1,:)=M(1,:);
V(2,:)=M(128,:);
dis=28.67/128;
di= t1;
for i=2:128
ti=di/28.67;% the ratio of distances, 
M(i,:)= [(1-ti)*V(1)+ti*V(2) (1-ti)*V(3)+ti*V(4)];
di= di+dis;
end

   
% assign elementwise optical coefficients - mus perturbation
ref = ones(nv,1)*refind;
mua = ones(nv,1)*mua_bkg;
mus = ones(nv,1)*mus_bkg;
prm.basis.hBasis = toastBasis(mesh,[bx by],'Linear_v2');
prm.bmua = prm.basis.hBasis.Map('M->B',mua);
prm.bmus = prm.basis.hBasis.Map('M->B',mus);

mat_hom=reshape(prm.bmua,bx,by);
mats_hom=reshape(prm.bmus,bx,by);

mat_hom= abs(mat_hom');
mats_hom= abs(mats_hom');
%insert lession assuming that there is surface that marks the inclusion
% for i=2:size(regno)
% blobs = find(regidx == regno(i));
% mua(blobs) = 0.17;%0.053;%
% mus(blobs) =1.5;%1.24;
% 
% end
for i=2:size(regno)
            blobs = find(regidx == regno(i));
            les_mua=0.017;% mualesion(randperm(numel(mualesion), 1));
            les_mus=1.1;%muslesion(mualesion==les_mua);
             for p=1:size(blobs)
                 m=idx(blobs(p),:);
                 mua(m)=les_mua;
                 %mus(m)=0;
                 mus(m)=les_mus;

             end
        end
figure;
%plot();
mesh.Display(mua, 'range',[0.001,0.17]); axis on; title('\mu_a target');
%subplot(2,2,2); mesh.Display(mus, 'range',[0.8,2]);     axis off; title('\mu_s target');
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');
for yy=1:2
mesh.SetQM(Q(yy,:),M);
qvec = (mesh.Qvec('Neumann','Gaussian',10));%10000*
mvec = 1000*(mesh.Mvec('Gaussian',10,refind));
% solve FEM linear system
smat = dotSysmat(mesh, mua, mus, ref,f);
data(:,yy) = (mvec.' * bicgstab(smat,qvec));
end
MESVEC=reshape(data,256,1);
size(MESVEC);
if(loop==0)
    results= MESVEC;
    loop=1;
    'here'
 else
    results= cat(2,results,MESVEC);
    size(results)
 end
% if k > 10
%         break
% end
end
csvwrite('/Users/HKhanene/Documents/MATLAB/back_up/aRC/RC_newlow.csv',full(results),0,0);

 