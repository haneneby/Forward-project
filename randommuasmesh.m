%%%%%%
%this script was created to test different mesh types (cts,rand, linea) and
% save their results in order to  compared their measure (Y) in
%compare_result.m file
%%%%%

clear all
close all
% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.0033 ; % background absorption [1/mm]
mus_bkg = 1.06;    % background scattering [1/mm];
itrmax = 100;   % CG iteration limit
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.5;      % scattering coef of medium index range uapper limit
w= .00069;        % wave length
freq = cm/w
% load mesh
%mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/gmsh_test/circle_blob.msh','gmsh');
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/parser_test/breast_dim2_3.msh','gmsh');
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region; 
regno = unique(regidx);
ref = ones(ne,1)*refind;

%  % % assign canstant elementwise optical coefficients - mus perturbation
mua = ones(ne,1)*mua_bkg;
mus = ones(ne,1)*mus_bkg;

% assign random mu perturbation
mua = (b-a).*rand(ne,1)+a;
mus = (b1-a1).*rand(ne,1)+a1;
mua_ref=mua;
mus_ref=mus;
% assign linear mu perturbation
%[mua, mus] =linearparam(ne,mua_bkg,mus_bkg);

% assign  to lesions - mu perturbation
for i=2:size(regno)
blobs = find(regidx == regno(i)); % assuming that there is surface that marks the inclusion
%indx=cat(2,indx,blobs');
mua(blobs) =0.17; %(b-a).* rand(1) +a%2.*mua_bkg;%(b-a).* rand(1) +a*2;
mus(blobs) =1.5;%(b1-a1).* rand(1) +a1%2.*mus_bkg;
end

mesh.Display
k=0;
figure; mesh.Display(mua);%,'range',[0.001,1]);
% set(gca,'XTick',[]) % Remove the ticks in the x axis!
% set(gca,'YTick',[]) % Remove the ticks in the y axis
% set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
% filename = sprintf('meshabsorption%d', k);
%saveas(gcf,filename,'png')
figure; mesh.Display(mus);%,'range',[0.6,1.7]);

% define source and detector locations
mm=1;
if mm<9 
    Q(1,:)=[6 -27];
else
    Q(1,:)=[0 0];%[6 -27];
end
Q(2,:)=[60 0];
%Q(4);
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

for yy=1:2
mesh.SetQM(Q(yy,:),M);
% mesh.SetQM(Q,M);
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');


% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2, refind);
size(mua);
size(mus);

% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,freq,'EL');
Phi = K\qvec;
Y(:,yy) = mvec.' * Phi;
Y;
% for reference, also solve the homogeneous problem
mus = mus_ref;%ones(ne,1)*mus_bkg;
mua =mua_ref;% ones(ne,1)*mua_bkg;
smat = dotSysmat(mesh, mua, mus, ref,freq, 'EL');
data_homog (:,yy)= (mvec.' * (smat\qvec));
end
toastWriteVector('Brest_homo.dat', data_homog);
toastWriteVector('Brest_hetero.dat', Y);
ddiff= Y-data_homog;
figure;
plot(real(ddiff));
MESVEC= sum(Y,2);

logY = log(MESVEC);
lnamp = real(logY);
phase = imag(logY);
% Display sinogram
figure
subplot(1,2,1); 
imagesc(lnamp);
title('logarithmic amplitude');
xlabel('source index q');
ylabel('detector index m');
axis equal tight;
colorbar
subplot(1,2,2); 
imagesc(phase);
title('phase');
xlabel('source index q');
ylabel('detector index m');
axis equal tight;
colorbar

% Write solver results to file
data = reshape(log(abs(Y')),[],1);
%toastWriteVector('Brest_homo_1.dat', Y);%shape1_LESION_hom_far


% Display boundary profile
figure
subplot(1,2,1); 
% title(' Amplitude');
% xlabel('Detector index q');
plot (real(data_homog(:,1)));
hold on
plot (real(data_homog(:,2)));
title(' Amplitude');
xlabel('Detector index q');

ylabel('Intensity');
legend('s1','s2')

subplot(1,2,2); 
plot (imag(data_homog(:,1)));
hold on
plot (imag(data_homog(:,2)));
ylabel('Intensity');
title(' Phase');
xlabel('Detector index q');
legend('s1','s2')

% %%%%%%% function%%%%%%%
function [f, p]= linearparam(ne,mua_bkg,mus_bkg)
vecta= ones(ne,1); 
vects= ones(ne,1); 

j=1;
limit=fix(ne/100)
step=100
for  i=1:step:ne
    if (i<=step*limit)
         vecta (i:i+99)=mua_bkg+0.0001* j;
         vects (i:i+99)=mus_bkg+0.01* j;

    else
         vecta (i:ne)=mua_bkg+0.0001* j;
         vects (i:ne)=mus_bkg+0.01* j;

    end
    vecta(i+1)
    j=j+1;
end
f= vecta;
p=vects;
end
