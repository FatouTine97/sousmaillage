clear;
cd '/home/fatou/Documents/sous_maillage/TP';

% Chargement des deux fichiers
M1 = load('carto_t1.dat');
M2 = load('carto_t.dat');

nbc = 1;         % Nombre de composantes (Ez par exemple)
nuc = 1;         % Composante à visualiser
nbcouleur = 256; % Nombre de couleurs
nbl = 250;       % Nombre de lignes par image

% Calcul du nombre de pas de temps
nbt1 = floor(size(M1,1) / (nbl * nbc));
nbt2 = floor(size(M2,1) / (nbl * nbc));
nbt = min(nbt1, nbt2); % synchronisation : on prend le plus petit

fprintf('Nombre de pas de temps détectés : %d (min des deux fichiers)\n', nbt);

offset = nbl * nbc;

figure(1);
clf;

for j = 1:nbt
    nut = j;

    % --- Extraction pour carto_t1.dat ---
    imin1 = (nuc - 1) * nbl + 1 + (nut - 1) * offset;
    imax1 = (nut - 1) * offset + nuc * nbl;
    data1 = M1(imin1:imax1, :);

    % --- Extraction pour carto_t.dat ---
    imin2 = (nuc - 1) * nbl + 1 + (nut - 1) * offset;
    imax2 = (nut - 1) * offset + nuc * nbl;
    data2 = M2(imin2:imax2, :);

    % --- Normalisation commune ---
    global_max = max([max(abs(data1), [], 'all'), max(abs(data2), [], 'all')]);
    if global_max == 0, global_max = 1; end

    img1 = abs(data1) / global_max;
    img2 = abs(data2) / global_max;

    % --- Concaténation et affichage ---
    img_combined = [img1, img2];

    imagesc(img_combined, [0 1]);
    axis image;
    colormap(jet(nbcouleur));
    colorbar;

    title(['Time step: ', num2str(j), '   |   Gauche: carto\_t1.dat   |   Droite: carto\_t.dat']);

    pause(1);
end
