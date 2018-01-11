clear all
close all
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/gmsh_test/breast8.msh','gmsh');
mesh.Display
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region;
regno = unique(regidx);
regno
% for i=2:regno
% blobel(:,i) = find(regidx == regno(i)); % assuming that the second surface marks the inclusion
% end
blobel_1= find(regidx == regno(2)); % assuming that the second surface marks the inclusion
blobel_2= find(regidx == regno(3)); % assuming that the second surface marks the inclusion

refind = 1.4; % refractive index
c0 = 0.3; % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
mus_bkg = 1; % background scattering [1/mm];
ref = ones(nv,1)*refind;
mua = ones(nv,1)*mua_bkg;
mus = ones(nv,1)*mus_bkg;

% for i=2:regno
%     mus(blobel(:,i)) = mus_bkg*2;
% end
mus(blobel_1) = mus_bkg*2;
mus(blobel_2) = mus_bkg*0.8;

% Display the target distributions
figure('position',[0,0,640,420]);
subplot(2,3,1); mesh.Display(mua, 'range',[0.005,0.025]);
axis off; title('\mu_a target');
subplot(2,3,2); mesh.Display(mus, 'range',[0.8,2.2]);
axis off; title('\mu_s target');
grd = [32,32];
basis = toastBasis(mesh, grd);
bmua = basis.Map('M->B', mua);
length(bmua)

