% clear previous data
clear;
cd '/home/fatou/Documents/sous_maillage'/TP/observation;


% Lire les données depuis le fichier texte
M11 = load(['Ez_t.dat']);
M12 = load(['Ez_r1.dat']);
M13 = load(['Ez_r2.dat'])
%M14 = load(['ez_s.dat'])
% Extraction des colonnes
t1 = M11(:,1); % Temps
t2 = M12(:,1); % Temps
t3 = M13(:,1); % Temps
nufig = 10;

for i=1:5
    % Création de la figure
    figure(nufig); clf;
    hold on;
    % Tracer l3s colonnes 2 et 3 en rouge et vert , t, M11(:,3), 'g', t, M11(:,4), 'y')t2, M12(:,3), 'g',, t1, M11(:,3), 'r '
    plot(t1, M11(:,1+i), 'Color','b');
    plot(t2, M12(:,1+i), 'Color','g'); %
    plot(t3, M13(:,1+i), 'Color','r');
    
    % Activer la grille
    grid on; % le maillage 
    
    % Ajouter des titres et légendes si besoin
    xlabel('t ');% axe des x 
    chaine = sprintf('Amplitude de Ez (V/m), Point %d',i);

    ylabel(chaine); % axe des y
    title('Évolution temporelle de Ez'); % titre du figure 
    legend('Ez-GM-PM', 'Ez-R1','Ez-R2');% legendre

     nufig = nufig +1;
end

return;

