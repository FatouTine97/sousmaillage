clear; clc;
close all
root = '/home/fatou/Documents/sous_maillage/TP/observation';
f1 = fullfile(root,'carto_t1.dat');   % Ez (CG)
f2 = fullfile(root,'carto_t.dat');    % ez_s (FG)
outv = fullfile(root,'carto_compare.mp4');

assert(isfile(f1) && isfile(f2), 'Fichiers introuvables.');

[F1, T1] = read_carto_blocks(f1);
[F2, T2] = read_carto_blocks(f2);

% Même nb de frames
nbt = min(numel(F1), numel(F2));
F1  = F1(1:nbt);  F2 = F2(1:nbt);
T   = [];  % temps commun si dispo
if ~isempty(T1) && ~isempty(T2)
    T = arrayfun(@(k) max([T1(min(k,end)), T2(min(k,end))]), 1:nbt);
end

% Échelle commune globale
gmax = 0;
for k = 1:nbt
    gmax = max([gmax, max(abs(F1{k}),[],'all'), max(abs(F2{k}),[],'all')]);
end
if gmax==0, gmax=1; end

% Writer
try
    vw = VideoWriter(outv,'MPEG-4');
catch
    vw = VideoWriter(strrep(outv,'.mp4','.avi'),'Motion JPEG AVI');
end
vw.FrameRate = 10; vw.Quality = 95; open(vw);

fig = figure('Color','w','Position',[50 50 1100 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
colormap(jet(256));

for k = 1:nbt
    nexttile(1);
    imagesc(abs(F1{k})/gmax, [0 1]); axis image off; colorbar;
    ttl1 = 'carto\_t1 (CG)'; if ~isempty(T), ttl1 = sprintf('%s | t=%.6g s', ttl1, T(k)); end
    title(ttl1);

    nexttile(2);
    imagesc(abs(F2{k})/gmax, [0 1]); axis image off; colorbar;
    ttl2 = 'carto\_t (FG)';  if ~isempty(T), ttl2 = sprintf('%s | t=%.6g s', ttl2, T(k)); end
    title(ttl2);

    drawnow;
    writeVideo(vw, getframe(fig));
     pause(1); 
end

close(vw);
fprintf('Vidéo écrite : %s\n', outv);

% -------- même utilitaire que ci-dessus --------
function [frames, times] = read_carto_blocks(fname)
    L = readlines(fname);
    frames = {}; times = [];
    buf = strings(0); tcur = NaN;
    for i = 1:numel(L)
        s = strtrim(L(i));
        if strlength(s)==0, flush(); continue; end
        if startsWith(s,"#")
            tok = regexp(s,'^#\s*t\s*=\s*([+\-0-9.eE]+)','tokens','once');
            if ~isempty(tok), tcur = str2double(tok{1}); end
            continue;
        end
        buf(end+1,1) = s; %#ok<AGROW>
    end
    flush();
    keep = ~cellfun(@isempty, frames); frames = frames(keep); times = times(keep);
    function flush()
        if ~isempty(buf)
            M = str2num(strjoin(buf, newline)); %#ok<ST2NM>
            frames{end+1} = M; times(end+1) = tcur;
            buf = strings(0); tcur = NaN;
        end
    end
end
