
clear;
cd '/home/fatou/Documents/sous_maillage/TP/observation';

M = load('carto_t1.dat');
nbc = 1;      % Nombre de composantes (ex: Ez)
nuc = 1;      % Composante à visualiser
nbcouleur = 256;
% Paramètres de la grille (doivent correspondre à ceux du Fortran)
nt_ech = 12
nx_ech = 300/2;
ny_ech = 300/2;

s = size(M);
nx_ech = s(1)/nt_ech;
ny_ech = s(2);
nbl = nx_ech;        % Nombre de lignes par image (Nx_sm / 2)
nbt = nt_ech;          % Compteur de pas de temps
taille = size(M);
total_lines = taille(1);

% Fixe la taille d'image désirée
%%nbt = floor(total_lines / (nbl * nbc));
fprintf('Nombre de pas de temps détectés : %d\n', nbt);

offset = nbl * nbc;

figure(1);
clf;
colormap(jet(nbcouleur));

for j = 1:nbt
    nut = j;
    %imin = (nuc - 1) * nbl + 1 + (nut - 1) * offset;
    %imax = (nut - 1) * offset + nuc * nbl;
    imin =  1 + (nut - 1) * offset;
    imax = nut * offset ;

    if imax > total_lines
        warning('Pas de temps %d dépasse les limites du fichier.', j);
        break;
    end

    M1 = M(imin:imax, :); % extrait sous-matrice

    maximum = max(abs(M1), [], 'all');
    if maximum == 0
        maximum = 1;
    end

    M2 = nbcouleur * abs(M1) / maximum;

    imagesc(M2);
    axis image;
    colorbar;
    title(['Time step: ', num2str(j)]);

    pause(1);
end



