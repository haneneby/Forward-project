function nonuniqueness

close all
clear all

% some parameters
refind = 1;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.002; % background absorption [1/mm]
mus_bkg = 1.06;    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion
w= .00069;        % wave length
f = cm/w
% load mesh% load mesh
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/parser_test/breast_dim2_3.msh','gmsh');
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region;
regno = unique(regidx);
%blobel = find(regidx == regno(2));
% assuming that there is surface that marks the inclusion

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
qvec = real(mesh.Qvec('Neumann','Gaussian',2));
mvec = real(mesh.Mvec('Gaussian',2,refind));

% assign elementwise optical coefficients - mus perturbation
ref = ones(ne,1)*refind;
mua = ones(ne,1)*mua_bkg;
mus = ones(ne,1)*mus_bkg;
%mus(blobel) = mus_bkg*2;
for i=2:size(regno)
blobs = find(regidx == regno(i));
mua(blobs) =0.017; 
%mus(blobs) =mus_bkg*2;
end
figure('position',[0,0,640,420]);
subplot(2,2,1); mesh.Display(mua, 'range',[0.005,0.25]); axis off; title('\mu_a target');
subplot(2,2,2); mesh.Display(mus, 'range',[0.8,2.2]);     axis off; title('\mu_s target');
hold on

plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');

% solve FEM linear system
smat = dotSysmat(mesh, mua, mus, ref, 'EL');
data(:,yy) = (mvec.' * (smat\qvec));
%f=450000;
% for reference, also solve the homogeneous problem
mus = ones(ne,1)*mus_bkg;
mua = ones(ne,1)*mua_bkg;
smat = dotSysmat(mesh, mua, mus, ref,f, 'EL');
data_homog (:,yy)= (mvec.' * (smat\qvec));
end
toastWriteVector('Brest_homo.dat', data_homog);
toastWriteVector('Brest_hetero.dat', data);

size (data)
size (data_homog)
ddiff= -data+data_homog;
figure;
plot(real(ddiff));%, [-0.26,0.015]); axis equal tight;
 hold on
plot(imag(ddiff));

title('target data');


