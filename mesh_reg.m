clear all
close all
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/dataset/dataset_mesh/260.msh','gmsh');
mesh.Display
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region;
regno = unique(regidx);
regno
% for i=2:regno
% blobel(:,i) = find(regidx == regno(i)); % assuming that the second surface marks the inclusion
% end
% blobel_1= find(regidx == regno(2)); % assuming that the second surface marks the inclusion
% blobel_2= find(regidx == regno(3)); % assuming that the second surface marks the inclusion
[vtx,idx,eltp] = mesh.Data


% blobel_1= find(regidx == regno(2))
% m=idx(blobel_1(1),:)
refind = 1.4; % refractive index
c0 = 0.3; % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.002; % background absorption [1/mm]
mus_bkg = 1; % background scattering [1/mm];
ref = ones(nv,1)*refind;
mua = ones(nv,1)*mua_bkg;
mus = ones(nv,1)*mus_bkg;
% mua(m)=30

for i=2:size(regno)
        blobs = find(regidx == regno(i));
        for j=1:size(blobs)
            m=idx(blobs(j),:)
            mua(m)=0.017
        end
end
% for i=2:regno
%     mus(blobel(:,i)) = mus_bkg*2;
% end
%   mua(blobel_1) =3*mua_bkg;
%  mua(blobel_2) = 15;

% Display the target distributions
% figure('position',[0,0,640,420]);
% subplot(2,3,1); mesh.Display(mua, 'range',[0.005,0.025]);
% axis off; title('\mu_a target');
% subplot(2,3,2); mesh.Display(mus, 'range',[0.8,2.2]);
% axis off; title('\mu_s target');
% grd = [10,10];
% basis = toastBasis(mesh, grd,'Linear_v2');
% bmua = basis.Map('M->B', mua);
% length(bmua)
% grd = [32,32];
bx=128; by=128;
prm.basis.hBasis = toastBasis(mesh,[bx by],'Linear_v2');

prm.bmua = prm.basis.hBasis.Map('M->B',mua)
mat=reshape(prm.bmua,bx,by)
mat= abs(mat')
h = imagesc((mat)); axis xy equal tight off % rot90
max_image = max(abs(mat(:)))
min_image = min(abs(mat(:)))

% prm.bmus = prm.basis.hBasis.Map('M->B',mus);
% axes(handles.axes1);
