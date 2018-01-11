clear all
close all
mesh = toastMesh('/Users/HKhanene/Documents/MATLAB/gmsh_test/circle_blob.msh')%,'gmsh');
%mesh.Display
% mesh.Write('mesh_mod.msh');
% mesh = toastMesh('circle_blob_mod.msh');
%mesh = toastMesh('breast_mod.msh');
mesh.Display
%rad = 10;

% 
% % Create the mesh
% rad = 25;
% nsect = 6;
% nring = 32;
% nbnd = 2;
% [vtx,idx,eltp] = mkcircle (rad, nsect, nring, nbnd);
% mesh = toastMesh(vtx,idx,eltp);
% 
% Load bitmaps for parameter distributions
bmua = imread('demo_matlab_fwd2_mua.png');
bmus = imread('demo_matlab_fwd2_mus.png');

% Scale to desired parameter ranges
bmua = double(bmua)./255.*0.02 + 0.01;
bmus = double(bmus)./255.*1.0 + 1.0;
% different way %
figure;
subplot(1,2,1);
imagesc(rot90(bmua));
axis equal tight; colorbar
title('\mu_a');
subplot(1,2,2);
imagesc(rot90(bmus));
axis equal tight; colorbar
title('\mu_s');

% Map to mesh basis
grd = size(bmua);
basis = toastBasis (mesh,grd);
mua = basis.Map ('B->M',bmua);
mus = basis.Map ('B->M',bmus);
figure; mesh.Display(mua);
figure; mesh.Display(mus);

nnd = mesh.NodeCount;
ref_bkg = 1.4;
ref = ones(nnd,1) * ref_bkg;
% Create the source and detector positions
nq = 2;
% for i=1:nq
%   phi_q = 2*pi*(i-0.5)/nq;
%   Q(i,:) = (0.2* rad*i-0.5* rad)  * [ sin(phi_q)-4 cos(phi_q) ];%[cos(phi_q) (2*cos(phi_q)-0.5* rad ) *sin(phi_q)];% -0.5* rad * [cos(phi_q)- 1 sin(phi_q)];
% end

% diagonal
Q(1,:)=[0 -5];
Q(2,:)=[10 0];

% horizontal
% Q(1,:)=[0 10];
% Q(2,:)=[10 10];
Q
nq = 128;
for i=1:nq
  % phi_q = 2*pi*(i-1)/nq;
  % Q(i,:) = -0.5* rad * [cos(phi_q)- 1 sin(phi_q)];
%   phi_m = 2*pi*(i-0.5)/nq;
%   M(i,:) = -0.5*rad * [cos(phi_m) 1]; %sin(phi_m)];
 M(i,:) = [ (i+35)/20 0.5*((i+35)/20)-5];% diagonal
 %M(i,:) = [ (i+35)/20 10];% horizontal
end
M
mesh.SetQM(Q,M);
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','m');


% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2, ref);



% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,0);
Phi = K\qvec;
Y = mvec.' * Phi;
Y
% Display sinogram
figure
imagesc(log(Y));
xlabel('source index q');
ylabel('detector index m');
axis equal tight;
colorbar


% % Display boundary profile
% figure
% hold on
% angle = [360/32:360/16:360];
% for i=1:size(Y,2)
%     ywrap = [Y(i:end,i); Y(1:i-1,i)];
%     plot(angle,log(ywrap),'o-');
% end
% axis([0 360 -14 -2]);
% xlabel('angular source-detector separation');
% ylabel('log intensity');

% Write solver results to file
data = reshape(log(Y'),[],1);
toastWriteVector('breast.dat', Y);

% % Show differences to homogeneous results
% data_homog = toastReadVector('demo_matlab_fwd1.dat');
% logYhomog = reshape(data_homog,nq,nq)';
% dlogY = log(Y)-logYhomog;
% 
% figure
% imagesc(dlogY);
% xlabel('source index q');
% ylabel('detector index m');
% axis equal tight;
% colorbar
% 
% figure
% hold on
% for i=1:size(Y,2)
%     ywrap = [dlogY(i:end,i); dlogY(1:i-1,i)];
%     plot(angle,ywrap,'o-');
% end
% axis([0 360 -1 0.1]);
% xlabel('angular source-detector separation');
% ylabel('log intensity perturbation');
