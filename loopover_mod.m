%%%%%% Training dataset creation script: new version after data analysis done in Jan 2018  %%%%%%5
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
w= [.00069 .00075 .0008 .00085];        % wave length
mualesion= [0.017 0.053 0.091 0.104 0.17 0.3]; %absorption coef of lesion
muslesion= [1.1  1.24 1.3 1.33 1.5 2];  % scattering coef of lesion
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.05;      % scattering coef of medium index range uapper limit
tol = 1e-6;     % tolerance of bicgstab where the default is  1e-6.
maxit= 20 ;     % max number of iteration, default is min (n,20) where n is the size of vector b in bicgstab(A,b)
%f = cm/w
count=1100 ;% figure index 3240
loop=0;
results= zeros([128 1]);
data= zeros([128 2]);
MESVEC= zeros([128 1]);
pt=1;
freq = zeros([1 4]);
for i = 1:length(w)
 freq(1,i) = cm/w(1,i);   
end

%% load mesh from directory
myDir = '/Users/HKhanene/Documents/MATLAB/dataset/balanced_data';%uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.msh')); %gets all txt files in struct
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  num = importdata(fullFileName);   %or readtable  
  mesh = toastMesh(fullFileName,'gmsh');
  %mesh.Display
  ll=0;
%% get mesh elements %
  ne = mesh.ElementCount;
  nv = mesh.NodeCount;
  regidx = mesh.Region;
  regno = unique(regidx);
%% define source_detector %
[Q , M]= source_detector();
%% assign elementwise optical coefficients - mus perturbation
ref = ones(ne,1)*refind;
mua = ones(ne,1)*mua_bkg;
mus = ones(ne,1)*mus_bkg;
for h = 1:1%length(freq)
    f= freq(h);
    
    for j=1:1%fix(length(mualesion)/2)
        
        % insert lession assuming that there is surface that marks the inclusion
        for i=2:size(regno)
        blobs = find(regidx == regno(i));
        les_mua=mualesion(randperm(numel(mualesion), 1));
        les_mus=muslesion(mualesion==les_mua);
        mua(blobs) =0.17;%les_mua;
        mus(blobs) =1.5;%les_mus;
        end
        %% Save absoption coeff figure for training
        %figure;
        figure('visible','off');
        mesh.Display(mua,'range',[0.002,0.3]); 
        set(gca,'XTick',[]) % Remove the ticks in the x axis!
        set(gca,'YTick',[]) % Remove the ticks in the y axis
        set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
        c = colorbar;
        c.Visible = 'off';
        %truesize(figure,[128 128]);
        filename = '/Users/HKhanene/Documents/MATLAB/dataset/balanced _training_data/abs%d-%d';
        filename = sprintf(filename,fix(f),count+j+ll);
       % saveas(gcf,filename,'jpg')
%         %% mus
%         figure; mesh.Display(mus,'range',[0.8,2]); 
%         axis off; title('\mu_s target');
%         hold on
%         plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
%         plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');

        %% solve forword model

        for yy=1:2
        mesh.SetQM(Q(yy,:),M);
        qvec = (mesh.Qvec('Neumann','Gaussian',2));
        mvec = (mesh.Mvec('Gaussian',2,refind));
        % solve FEM linear system
        smat = dotSysmat(mesh, mua, mus, ref,f, 'EL');
        data(:,yy) = (mvec.' * bicgstab(smat,qvec,tolCG,50));
        end
        %MESVEC=reshape(data,256,1);
        MESVEC= sum(data,2);

        if(loop==0)
            results= MESVEC;
            loop=1;
            'here'
         else
            results= cat(2,results,MESVEC);
            size(results);
        end
        
    end
    ll=ll+1;%fix(length(mualesion)/2);
end
count=count+1;%fix(length(mualesion)/2).*length(freq);

% save data measurments more frequently
if count== pt*120
    csvwrite('/Users/HKhanene/Documents/MATLAB/dataset/balanced _training_data_suite.csv',full(results),0,0);
    pt= pt+1;
end
end 
csvwrite('/Users/HKhanene/Documents/MATLAB/dataset/balanced _training_data_suite.csv',full(results),0,0);




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