function nonuniqueness

close all
clear all

% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.0033; % background absorption [1/mm]
mus_bkg = 1.06;    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion
w= .00069;        % wave length
a= 0.002;       % absorption coef of medium index range lower limit
b= 0.012;       % absorption coef of medium index range upper limit
a1= 0.6 ;      % scattering coef of medium index range lower limit
b1= 1.05;      % scattering coef of medium index range uapper limit
tol = 1e-6     % tolerance of bicgstab where the default is  1e-6.
maxit= 20      % max number of iteration, default is min (n,20) where n is the size of vector b in bicgstab(A,b)
f = cm/w
% load mesh% load mesh
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/parser_test/bRC/breast_dim2_10_0.msh','gmsh');
display(mesh)
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region;
regno = unique(regidx);

% define source 
mm=10;
if mm<9 
    Q(1,:)=[6 -27];
else
    Q(1,:)=[0 0];%[6 -27];
end
Q(2,:)=[60 0];
%define detector locations
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
qvec = (mesh.Qvec('Neumann','Gaussian',2));
mvec = (mesh.Mvec('Gaussian',2,refind));

% assign elementwise optical coefficients - mus perturbation
ref = ones(ne,1)*refind;
mua = ones(ne,1)*mua_bkg;
mus = ones(ne,1)*mus_bkg;
% assign random mu perturbation
% mua = (b-a).*rand(ne,1)+a;
% mus = (b1-a1).*rand(ne,1)+a1;
mua_ref=mua;
mus_ref=mus;



% % for reference, also solve the homogeneous problem
smat = dotSysmat(mesh, mua, mus, ref,f, 'EL');
data_homog (:,yy)= (mvec.' * bicgstab(smat,qvec,tol));

toastWriteVector('Brest_homo.dat', data_homog);

%insert lession assuming that there is surface that marks the inclusion
for i=2:size(regno)
blobs = find(regidx == regno(i));
mua(blobs) =0.7; 
mus(blobs) =1.5;
end
figure('position',[0,0,640,420]);
subplot(2,2,1); mesh.Display(mua, 'range',[0.001,0.17]); axis off; title('\mu_a target');
subplot(2,2,2); mesh.Display(mus, 'range',[0.8,2]);     axis off; title('\mu_s target');
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');

% solve FEM linear system
smat = dotSysmat(mesh, mua, mus, ref,f, 'EL');
data(:,yy) = (mvec.' * bicgstab(smat,qvec,tol));
end


toastWriteVector('Brest_hetero.dat', data);

% ddiff= data-data_homog;
figure;
plot(real (data_homog));
figure;
plot(real (data));
% hold on
% plot(imag (data));
%  hold on
% plot(imag(ddiff));

title('target data');


