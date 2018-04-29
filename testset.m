%%%%%% Training dataset creation script: new version after data analysis done in Jan 2018  %%%%%%5
close all
clear all

%% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.0022; % background absorption [1/mm]
mus_bkg = 1.1;    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion
w= [.00069 .00075 .0008 .00085];        % wave length
mualesion= [ 0.1  0.13  0.17  0.19] %[ 0.06 0.08 0.104  ]; %absorption coef of lesion
muslesion= [ 2.2  2.3  2.4  2.5 ]%[ 1.9 2 2.15  ];  % scattering coef of lesion
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.05;      % scattering coef of medium index range uapper limit
tol = 1e-6;     % tolerance of bicgstab where the default is  1e-6.
maxit= 20 ;     % max number of iteration, default is min (n,20) where n is the size of vector b in bicgstab(A,b)
%f = cm/w
count=0 ;% figure index 3240
firstloop=0;
data= zeros([128 2]);
data_hom= zeros([128 2]);

MESVEC= zeros([256 1]);
RC= zeros([256 1]);
hom_MESVEC=zeros([256 1]);

Relative= zeros([256 1]);
results= zeros([256 1]);
flfr= [690,750,800,850]
pt=1;
freq = zeros([1 4]);
for i = 1:length(w)
 freq(1,i) = cm/w(1,i);   
end
for fl=2:2% length(flfr)
  
filename = '/Users/HKhanene/Documents/ML_dataset/phantomData/%d/maps/absmat/abs%d-%d.csv';%F800/absmat/abs%d-%d.csv';
filename1 = '/Users/HKhanene/Documents/ML_dataset/phantomData/%d/maps/scamat/sca%d-%d.csv';%F800/scamat/sca%d-%d.csv';
bx=128; by=128;


%% load mesh from directory
path='/Users/HKhanene/Documents/ML_dataset/phantomData/probmesh';%/shapes';%uigetdir; %gets directory

path = sprintf(path,flfr(fl));
myDir =path;%/shapes';%uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.msh')); %gets all txt files in struct
myFiles = natsortfiles({myFiles.name});

for k = 1:length(myFiles)%800%
  baseFileName = myFiles(k);%.name;
  fullFileName = fullfile(myDir,'/', baseFileName);
  fprintf(1, 'Now reading %s\n', char(fullFileName));
%   num = importdata(fullFileName);   %or readtable  
  mesh = toastMesh(char(fullFileName),'gmsh');
  %mesh.Display
  
  ll=0;
%% get mesh elements %
  ne = mesh.ElementCount;
  nv = mesh.NodeCount;
  regidx = mesh.Region;
  regno = unique(regidx);
  [vtx,idx,eltp] = mesh.Data;

%% define source_detector %
%     [Q , M]= source_detector();

%% assign elementwise optical coefficients 
    ref = ones(nv,1)*refind;
    mua = ones(nv,1)*mua_bkg;
    mus = ones(nv,1)*mus_bkg;


    
    f= freq(fl);
    if (k==1)
        if (fl==1 )
             j=0;
        else
          j=j+length(myFiles);   

        end
   
    end
    
    % insert lession assuming that there is surface that marks the inclusion
    for i=2:size(regno)
        blobs = find(regidx == regno(i));
        les_mua= 0.126;%mualesion(j);% mualesion(randperm(numel(mualesion), 1));
        les_mus= 1.1;%muslesion(j);%muslesion(mualesion==les_mua);
         for p=1:size(blobs)
             m=idx(blobs(p),:);
             mua(m)=les_mua-mua_bkg;
             %mus(m)=0;
             mus(m)=les_mus-mus_bkg;

         end
    end
 %% Save absoption coeff matrix for MLtraining
        prm.basis.hBasis = toastBasis(mesh,[bx by],'LINEAR_V2');
        prm.bmua = prm.basis.hBasis.Map('M->B',mua);
        prm.bmus = prm.basis.hBasis.Map('M->B',mus);

        mat=reshape(prm.bmua,bx,by);
        mats=reshape(prm.bmus,bx,by);
        
        mat= abs(mat');
        mats= abs(mats');
        %imagesc (mat)
    
    %imagesc (mat_hom)
    file = sprintf(filename,flfr(fl),flfr(fl),k+j);
    file1 = sprintf(filename1,flfr(fl),flfr(fl),k+j);
%         imwrite((mat_hom),file,'mode', 'lossless');
%         imwrite((mats_hom),file1,'mode', 'lossless');
        %imwrite(YourImage, 'YourFile.jpg', 'mode', 'lossless')
        
%% save maps        
%     csvwrite(file,full(mat),0,0);
%     csvwrite(file1,full(mats),0,0);
% %   
%% define source_detector %
    [Q , M]= source_detector();
 %% solve the dot.forword model (projection)
        for yy=1:2
        mesh.SetQM(Q(yy,:),M);
        qvec = (mesh.Qvec('Isotropic','Gaussian',2));%  Neumann or Isotropic for boundry condition
        mvec = 1000*(mesh.Mvec('Gaussian',2,refind));
        % solve FEM linear system
        smat = dotSysmat(mesh, mua, mus, ref,f);
        data(:,yy) = (mvec.' * bicgstab(smat,qvec,tolCG,50));
        end
        MESVEC=reshape(data,256,1);
        
        % save results
         if(k==1)
             results= MESVEC;
         else
             results= cat(2,results,MESVEC);
         end
end 
filel='/Users/HKhanene/Documents/ML_dataset/phantomData/%d/maps/syn_measures%d.csv';%/shapes';%uigetdir; %gets directory

filel = sprintf(filel,flfr(fl),flfr(fl));
csvwrite(filel,full(results),0,0);
end

function [Q , M]= source_detector()
% define source %
  Q(1,:)=[0 0];%[6 -27];
  Q(2,:)=[60 0];
% define detector %
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
end