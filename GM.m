
clear;
cd '/home/fatou/Documents/sous_maillage/TP';
M = load('carto_t1.dat');

nbc = 1;      % Nombre de composantes (ex: Ez)
nuc = 1;      % Composante à visualiser
nbcouleur = 256;

taille = size(M);
total_lines = taille(1);

% Fixe la taille d'image désirée
nbl = 250;

% Calcule le nombre de pas de temps maximum possibles
nbt = floor(total_lines / (nbl * nbc));
fprintf('Nombre de pas de temps détectés : %d\n', nbt);  

offset = nbl * nbc;

figure(1);
clf;
colormap(jet(nbcouleur));

for j = 1:nbt
    nut = j;
    imin = (nuc - 1) * nbl + 1 + (nut - 1) * offset;
    imax = (nut - 1) * offset + nuc * nbl;

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
