clear all
close all
rad = 25;   % mesh radius [mm]
nsect = 6;  % number of sectors
nring = 32; % number of rings
nbnd = 2;   % number of boundary rings

[vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd);
            % create the mesh geometry
mesh = toastMesh (vtx,idx,eltp);
            % create the mesh object
            [vtx,idx,eltp] = mesh.Data;
nnode = mesh.NodeCount;
nel = mesh.ElementCount;
mesh.Display
mua_bkg = 0.01;
mus_bkg = 1.0;
ref_bkg = 1.4;
nnd = mesh.NodeCount;
mua = ones(nnd,1) * mua_bkg;
mus = ones(nnd,1) * mus_bkg;
ref = ones(nnd,1) * ref_bkg;
nq = 16;
for i=1:nq
  phi_q = 2*pi*(i-1)/nq;
  Q(i,:) = rad * [cos(phi_q) sin(phi_q)];
  phi_m = 2*pi*(i-0.5)/nq;
  M(i,:) = rad * [cos(phi_m) sin(phi_m)];
end
mesh.SetQM(Q,M);
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bx','MarkerFaceColor','b');
qvec = mesh.Qvec('Neumann','Gaussian',2);
%mvec = mesh.Mvec('Gaussian',2);