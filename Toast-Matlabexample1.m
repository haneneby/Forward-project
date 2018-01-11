rad = 25;   % mesh radius [mm]
nsect = 6;  % number of sectors
nring = 32; % number of rings
nbnd = 2;   % number of boundary rings

[vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd);
            % create the mesh geometry