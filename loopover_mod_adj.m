%%%%%% Training dataset creation script: new version after data analysis done in Jan 2018  %%%%%%5
close all
clear all

%% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = [0.002 0.0022 0.0019  0.0042]; % background absorption [1/mm]
mus_bkg = [1.5 1.2 0.9  0.82 ];    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion
w= [.00069 .00075 .0008 .00085];        % wave length
mualesion= [0.125  0.13  0.135  0.14 0.145  0.155 0.16]; %[0.12 0.134 ] ;%absorption coef of lesion  [ 0.17 0.126 0.014 0.017  ]  [ 0.121 0.123 0.125 0.127 0.129  0.131  0.133  ];
muslesion= [ 1 1.02 1.04 1.06 1.08 1 1.1 ];%[0.93 0.85  ];%[1 1.3];%1.1; % scattering coef of lesion [ 1  1.1 1.12 1.15 1.17 1.2 1.23  ];
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.05;      % scattering coef of medium index range uapper limit
tol = 1e-6;     % tolerance of bicgstab where the default is  1e-6.
maxit= 20 ;     % max number of iteration, default is min (n,20) where n is the size of vector b in bicgstab(A,b)
%f = cm/w
count=0 ;% figure index 5480
firstloop=0;
data= zeros([128 2]);
data_hom= zeros([128 2]);

MESVEC= zeros([256 1]);
RC= zeros([256 1]);
hom_MESVEC=zeros([256 1]);

Relative= zeros([256 1]);
results= zeros([256 1]);

pt=1;
freq = zeros([1 4]);
for i = 1:length(w)
 freq(1,i) = cm/w(1,i);   
end
filename = '/Users/HKhanene/Documents/ML_dataset/trainset750/absmat/abs%d-%d.csv';
filename1 = '/Users/HKhanene/Documents/ML_dataset/trainset750/scamat/sca%d-%d.csv';
bx=128; by=128;


%% load mesh from directory
myDir = '/Users/HKhanene/Documents/ML_dataset/phantomData/lightmesh';%uigetdir; %gets directory
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
    [Q , M]= source_detector();


for h = 2:2%length(freq)
    
    f= freq(h);
    %% assign elementwise optical coefficients 
    ref = ones(nv,1)*refind;
    mua = ones(nv,1)*mua_bkg(h);
    mus = ones(nv,1)*mus_bkg(h);
  
    
    %% insert lesion and get projection measures
    for j=1:(length(mualesion))
        
        % insert lession assuming that there is surface that marks the inclusion
        for i=2:size(regno)
            blobs = find(regidx == regno(i));
            les_mua=mualesion(j);% (mualesion (2)-mualesion(1)).*rand(1) + mualesion(1);% mualesion(randperm(numel(mualesion), 1));
            les_mus=muslesion(j);%(muslesion (2)-muslesion(1)).*rand(1) + muslesion(1);%muslesion (j);% (muslesion (2)-muslesion(1)).*rand(1) + muslesion(1);% muslesion;%muslesion(mualesion==les_mua);
             for p=1:size(blobs)
                 m=idx(blobs(p),:);
                 mua(m)=les_mua-mua_bkg(h);
                 %mus(m)=0;
                 mus(m)=les_mus-mus_bkg(h);

             end
        end
        %% Save absoption coeff matrix for MLtraining
 %% map homogenious mesh to grid and save coef matrix for ML training
        prm.basis.hBasis = toastBasis(mesh,[bx by],'LINEAR_V2');
        
        prm.bmua = prm.basis.hBasis.Map('M->B',mua);
        prm.bmus = prm.basis.hBasis.Map('M->B',mus);

        mat=reshape(prm.bmua,bx,by);
        mats=reshape(prm.bmus,bx,by);
        
        mat= rot90 (abs(mat));
        mats=rot90(abs(mats));
        %imagesc (mat)

        file = sprintf(filename,fix(f),count+j+ll);
        file1 = sprintf(filename1,fix(f),count+j+ll);

%         figure;
%         hoo = imagesc((mat)); axis xy equal tight off % rot90
%         figure;
%         hoo = imagesc((mats)); axis xy equal tight off % rot90
%        
        csvwrite(file,full(mat),0,0);
        csvwrite(file1,full(mats),0,0);
        
%         imwrite((mat),file);
%         imwrite((mats),file1);

%         %% plot source and detector position
%         figure; mesh.Display(mus,'range',[0.8,2]); 
%         axis off; title('\mu_s target');
%         hold on
%         plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
%         plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');

        %% solve the dot.forword model (projection)
        for yy=1:2
        mesh.SetQM(Q(yy,:),M);
        qvec = (mesh.Qvec('Isotropic','Gaussian',2));
        mvec = 1000*mesh.Mvec('Gaussian',2,refind);
        % solve FEM linear system
        smat = dotSysmat(mesh, mua, mus, ref,f);
        data(:,yy) = (mvec.' * bicgstab(smat,qvec,tolCG,50));
        end
        MESVEC=reshape(data,256,1);
        %MESVEC= sum(data,2);
        
%         % calculate the RC vector
%         RC = (hom_MESVEC-MESVEC)./hom_MESVEC;
%         RCvect=reshape(RC,256,1);

        % save results
        if(firstloop==0)
            results= MESVEC;
            firstloop=1;
            'here'
        else
        results= cat(2,results,MESVEC);
%         relative= cat(2,relative,RC);
        end
        
    end
    ll=ll+1%fix(length(mualesion));
end
count=count+length(mualesion);%length(freq)

% save data measurments more frequently
if count== pt*70
    csvwrite('/Users/HKhanene/Documents/ML_dataset/trainset750/measures.csv',full(results),0,0);
%     csvwrite('/Users/HKhanene/Documents/ML_dataset/phan2/Relativechange.csv',full(relative),0,0);

    pt= pt+1;
end
end 
csvwrite('/Users/HKhanene/Documents/ML_dataset/trainset750/measures.csv',full(results),0,0);
% csvwrite('/Users/HKhanene/Documents/ML_dataset/F800/Relativechange.csv',full(relative),0,0);



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