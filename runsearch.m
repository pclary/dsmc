function [findtimes, targets, phi2, mur, c, xt, yt] = runsearch(settings, ...
    outputsettings, ax, aborthandle, pausebutton)
%RUNSEARCH Runs a search using the specified algorithm
%   Returns the time at which each target was found and a vector of the
%   number of remaining targets at each detection time
%   
%   Only the settings struct is required

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 12/18/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack settings
xlim = settings.xlim;
ylim = settings.ylim;
vx = settings.vx;
vy = settings.vy;
h = settings.h;
mu = settings.mu;
nsamplepts = settings.nsamplepts;
mures = settings.mures;
cres = settings.cres;
umax = settings.umax;
algorithm = settings.algorithm;
ntargets = settings.ntargets;
findradius = settings.findradius;
findtimeconst = settings.findtimeconst;
maxtime = settings.tstop;
lambda = settings.lambda;
stopallfound = settings.stopallfound;
au = settings.agentuncertainty;
tu = settings.targetuncertainty;
spherical = settings.spherical;
starttime = settings.starttime;
datetitle = settings.datetitle;
xti = settings.xti;
yti = settings.yti;
xland = settings.xland;
yland = settings.yland;

if nargin < 2
    outputsettings = [];
end
if nargin < 3
    ax = struct();
end
if nargin < 4
    aborthandle = [];
end
if nargin < 5
    pausebutton = [];
end

% Set up output
if ~isempty(outputsettings)
    outfig = figure;
    axoutput = gca;
    set(outfig, 'Renderer', 'zbuffer');
    set(outfig, 'Visible', 'off');
    
    % Append (1), (2), etc to output file names to avoid overwriting old output
    if ~isempty(outputsettings)
        outputsettings.main.filename = processfilename(outputsettings.main.filename, outputsettings.main.animation);
        outputsettings.convergence.filename = processfilename(outputsettings.convergence.filename, outputsettings.convergence.animation);
        outputsettings.mu.filename = processfilename(outputsettings.mu.filename, outputsettings.mu.animation);
        outputsettings.coverage.filename = processfilename(outputsettings.coverage.filename, outputsettings.coverage.animation);
        outputsettings.mesohyperbolicity.filename = processfilename(outputsettings.mesohyperbolicity.filename, outputsettings.mesohyperbolicity.animation);
        
        if ~outputsettings.overwrite
            outputsettings.main.filename = getunusedfilename(outputsettings.main.filename);
            outputsettings.convergence.filename = getunusedfilename(outputsettings.convergence.filename);
            outputsettings.mu.filename = getunusedfilename(outputsettings.mu.filename);
            outputsettings.coverage.filename = getunusedfilename(outputsettings.coverage.filename);
            outputsettings.mesohyperbolicity.filename = getunusedfilename(outputsettings.mesohyperbolicity.filename);
        end
    end
end

% Sample distribution to get particles and targets
[mux, muy] = sampledist(mu, nsamplepts, xlim, ylim);
[muxtar, muytar] = sampledist(mu, ntargets, xlim, ylim);
foundtargets = [];
findtimes = 0;
findchance = 1 - exp(-h/findtimeconst);
phi2 = [];

% Set initial positions
xt = xti;
yt = yti;

% Initialize K arrays and coverage coefficients
[Kx, Ky] = meshgrid(0:cres-1, 0:cres-1);
Las = 1./(1 + Kx.^2 + Ky.^2).^(3/2);
cks = pts2dct(xt, yt, cres, xlim, ylim);

% Misc setup
fcount = [0, 0, 0, 0, 0];
stepnum = size(xt, 1);
t = 0;
co = get(gca,'ColorOrder');
clear lawnmowerstep

% Used for abort button signalling with the GUI
if ~isempty(aborthandle)
    abort = @() getappdata(aborthandle, 'Abort');
else
    abort = @() 0;
end

% Used for pause button signalling with the GUI
if ~isempty(pausebutton)
    paused = @() get(pausebutton, 'Value');
else
    paused = @() 0;
end

% Step through simulation until a stop condition is met
while t < maxtime && (~stopallfound || numel(foundtargets) < ntargets) && ~abort()
    
    % Step particles, targets, and trajectories forward in time according 
    % to the velocity field
    [mux, muy] = rk4step(t, h, mux, muy, vx, vy);
    [muxtarn, muytarn] = rk4step(t, h, muxtar, muytar, vx, vy);
    tdisps = sqrt((muxtarn-muxtar).^2 + (muytarn-muytar).^2);
    muxtar = muxtarn + tu*rand(size(muxtar)).*tdisps;
    muytar = muytarn + tu*rand(size(muxtar)).*tdisps;
    [xt, yt] = rk4step(t, h, xt, yt, vx, vy);
    
    muks0 = pts2dct(mux, muy, mures, xlim, ylim);
    muks = logmuks(muks0, lambda, cres, xlim, ylim);
    
    % Use search algorithm to determine new agent positions
    xtn0 = xt(stepnum, :);
    ytn0 = yt(stepnum, :);
    if strcmp(algorithm, 'DSMC')
        [xtn, ytn] = dsmcstep(xtn0, ytn0, h, umax, cks, muks, cres, xlim, ylim, au, spherical);
    elseif strcmp(algorithm, 'Lawnmower')
        [xtn, ytn] = lawnmowerstep(xtn0, ytn0, h, umax, xlim, ylim, 3, spherical);
    elseif strcmp(algorithm, 'Random Walk')
        [xtn, ytn] = randomwalkstep(xtn0, ytn0, h, umax, xlim, ylim, spherical);
    end
    xt(stepnum+1, :) = xtn;
    yt(stepnum+1, :) = ytn;
    stepnum = stepnum + 1;
    
    cks = pts2dct(xt, yt, cres, xlim, ylim);
    
    t = t + h;
    
    % Keep track of discovered targets
    notfound = ~ismember(1:ntargets, foundtargets);
    dists = min(sqrt(bsxfun(@minus, xt(stepnum, :), muxtar).^2 + ...
        bsxfun(@minus, yt(stepnum, :), muytar).^2), [], 2);
    if numel(dists) == 0
        dists = zeros(0, 1);
    end
    found = find(dists < findradius & notfound' & rand(ntargets, 1) < findchance);
    foundtargets = [foundtargets; found];
    findtimes = [findtimes; t*ones(numel(found), 1)];
    
    % Calculate convergence metric
    sks = cks - muks;
    phi2n = sum(sum(Las.*sks.^2));
    phi2 = [phi2; phi2n];
    
    % Plot particles and trajectories
    if isfield(ax, 'main') && ~isempty(ax.main)
        plotmain(ax.main, mux, muy, muxtar, muytar, foundtargets, xt, yt, ...
            stepnum, ntargets, findradius, xland, yland, xlim, ylim, co);
        maketitle(ax.main, t, datetitle, starttime);
    end
    
    % Plot convergence metric
    if isfield(ax, 'convergence') && ~isempty(ax.convergence)
        plotconvergence(ax.convergence, phi2);
    end
    
    % Plot log density of particles
    if isfield(ax, 'mu') && ~isempty(ax.mu)
        plotmu(ax.mu, muks, xlim, ylim, cres);
    end
    
    % Plot trajectory density
    if isfield(ax, 'coverage') &&  ~isempty(ax.coverage)
        plotcoverage(ax.coverage, cks, xlim, ylim, cres);
    end
    
    drawnow;
    
    % Save output figures
    if ~isempty(outputsettings)
        plotfun = {...
            @(ax) plotmain(ax, mux, muy, muxtar, muytar, foundtargets, xt, yt, ...
            stepnum, ntargets, findradius, xland, yland, xlim, ylim, co), ...
            @(ax) plotconvergence(ax, phi2), ...
            @(ax) plotmu(ax, muks, xlim, ylim, cres), ...
            @(ax) plotcoverage(ax, cks, xlim, ylim, cres), ...
            @(ax) plotmesohyperbolicity(ax, vx, vy, xlim, ylim, t, outputsettings.mesohyperbolicity.T, cres)};
        
        osfigs = {outputsettings.main, outputsettings.convergence, ...
            outputsettings.mu, outputsettings.coverage, ...
            outputsettings.mesohyperbolicity};
        
        for i = 1:numel(osfigs)
            os = osfigs{i};
            if os.enable && t >= fcount(i)*os.rate
                set(outfig, 'Position', [0, 0, os.width, os.height]);
                func = plotfun{i};
                func(axoutput);
                maketitle(axoutput, t, datetitle, starttime);
                set(outfig, 'Visible', 'off');
                drawnow;
                if os.animation
                    frame = getframe(outfig);
                    im = frame2im(frame);
                    [imind, map] = rgb2ind(im,256);
                    if fcount(i) == 0
                        imwrite(imind, map, os.filename, 'DelayTime', 1/30, 'LoopCount', inf);
                    else
                        imwrite(imind, map, os.filename, 'DelayTime', 1/30, 'WriteMode', 'append');
                    end
                else
                    if ~exist(os.filename, 'file')
                        mkdir(os.filename);
                    end
                    frame = getframe(outfig);
                    im = frame2im(frame);
                    imwrite(im, [os.filename, '/', sprintf('%.3d.png', fcount(i) + 1)], 'png');
                end
                fcount(i) = fcount(i) + 1;
            end
        end
    end
    
    % Pause
    while paused()
        pause(0.05);
    end
end

set(gca, 'ColorOrder', co);

targets = ntargets:-1:ntargets-numel(findtimes)+1;

mur = idct2(muks);

c = idct2(cks);


function plotmain(ax, mux, muy, muxtar, muytar, foundtargets, xt, yt, ...
    stepnum, ntargets, findradius, xland, yland, xlim, ylim, co)

plot(ax, mux, muy, '.b', 'MarkerSize', 2);
set(ax, 'ColorOrder', co(2:end, :), 'NextPlot', 'replacechildren');
hold(ax, 'on');
plot(ax, muxtar, muytar, 'co', 'MarkerSize', 5, 'LineWidth', 2);
plot(ax, muxtar(foundtargets), muytar(foundtargets), 'mo', ...
    'MarkerSize', 5, 'LineWidth', 2);

plot(ax, xland, yland, 'k.');

plot(ax, [xt(1, :); xt(1:stepnum, :)], [yt(1, :); yt(1:stepnum, :)]);
plot(ax, [xt(stepnum, :); xt(stepnum, :)], ...
    [yt(stepnum, :); yt(stepnum, :)], '.', 'MarkerSize', 15);

if ntargets > 0
    theta = linspace(0, 2*pi)';
    circx = bsxfun(@plus, xt(stepnum, :), findradius*cos(theta));
    circy = bsxfun(@plus, yt(stepnum, :), findradius*sin(theta));
    plot(ax, circx, circy);
end

hold(ax, 'off');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);


function plotconvergence(ax, phi2)

loglog(ax, phi2);
title(ax, 'Convergence');
xlabel(ax, 'Steps');
ylabel(ax, '\phi^2');


function plotmu(ax, muks, xlim, ylim, cres)

mur = idct2(muks);
x = linspace(xlim(1), xlim(2), cres);
y = linspace(ylim(1), ylim(2), cres);
[x, y] = meshgrid(x, y);
surf(ax, x, y, mur, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Log(Particle distribution)');
caxis(ax, sort([0, max(max(mur))]));


function plotcoverage(ax, cks, xlim, ylim, cres)

c = idct2(cks);
x = linspace(xlim(1), xlim(2), cres);
y = linspace(ylim(1), ylim(2), cres);
[x, y] = meshgrid(x, y);
surf(ax, x, y, c, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Coverage density');
maxc = max(max(c));
if isnan(maxc)
    maxc = 1;
end
caxis(ax, [0, maxc]);


function plotmesohyperbolicity(ax, vx, vy, xlim, ylim, t, T, cres)

[mh, x0, y0] = mesohyperbolicity(vx, vy, xlim, ylim, t, T, cres);
surf(ax, x0, y0, mh, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Mesohyperbolicity');
colormap(ax, mhcolormap(256));
caxis(ax, [-2/(T)^2, 6/(T)^2]);
colorbar('peer', ax, 'EastOutside');


function plots(ax, Lasks, cres, xlim, ylim)

s = idct2(Lasks);
x = linspace(xlim(1), xlim(2), cres);
y = linspace(ylim(1), ylim(2), cres);
[x, y] = meshgrid(x, y);
surf(ax, x, y, s, 'EdgeColor', 'none');
title(ax, 'DSMC surface');


function plotgrad(ax, Lasks, cres, xlim, ylim)

x = linspace(xlim(1), xlim(2), 22);
x = x(2:end-1);
y = linspace(ylim(1), ylim(2), 22);
y = y(2:end-1);
[xa, ya] = meshgrid(x, y);

Bdir = @(Bs) -bsxfun(@rdivide, Bs, sqrt(sum((Bs).^2, 1)));

tmp = Bdir(B(xa(:), ya(:), Lasks, cres, xlim, ylim));
u = tmp(1, :);
v = tmp(2, :);

quiver(ax, xa(:)', ya(:)', u, v);
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'DSMC surface gradient');


function maketitle(ax, t, datetitle, starttime)

if datetitle
    title(ax, datestr(starttime + t/60/60/24, 0));
else
    title(ax, sprintf('t = %.2f s', t));
end


function filename = processfilename(filename, animation)

if animation && (numel(filename < 4) || ~strcmpi(filename(end-3:end), '.gif'))
    filename = [filename, '.gif'];
end


function filename = getunusedfilename(filename)
% Appends (1), (2), etc to filenames if the file already exists

originalfn = filename;

if numel(filename) > 0 && exist(filename, 'file')
    i = 1;
    while exist(filename, 'file')
        filename = originalfn;
        k = strfind(filename, '.');
        if numel(k) == 0
            k = numel(filename) + 1;
        end
        k = k(end);
        if k > 1
            filename = sprintf('%s (%d)%s', filename(1:k-1), i, filename(k:end));
        else
            filename = sprintf('(%d)%s', i, filename(k:end));
        end
        i = i + 1;
    end
end

