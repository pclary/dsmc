function [findtimes, targets, phi2, mu, c, xt, yt] = runsearch(settings, ...
    outputsettings, ax, aborthandle, pausebutton)
%RUNSEARCH Runs a search using the specified algorithm
%   Returns the time at which each target was found and a vector of the
%   number of remaining targets at each detection time
%   
%   Only the settings struct is required

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 3/10/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Unpack settings
xlim = settings.xlim;
ylim = settings.ylim;
vx = settings.vx;
vy = settings.vy;
h = settings.h;
mu = settings.mu;
nsamplepts = settings.nsamplepts;
substeps = settings.substeps;
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
regime = settings.regime;
xland = settings.xland;
yland = settings.yland;
verbose = settings.verbose;
showgrad = true;

if verbose
    writelog = @(s) fprintf('[%.2f] %s\n', toc(), s);
else
    writelog = @(s) 0;
end

writelog('Preparing to start...');

if nargin < 2
    outputsettings = [];
end
if nargin < 3
    ax = [];
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
        outputsettings.gradient.filename = processfilename(outputsettings.gradient.filename, outputsettings.gradient.animation);
        
        if ~outputsettings.overwrite
            outputsettings.main.filename = getunusedfilename(outputsettings.main.filename);
            outputsettings.convergence.filename = getunusedfilename(outputsettings.convergence.filename);
            outputsettings.mu.filename = getunusedfilename(outputsettings.mu.filename);
            outputsettings.coverage.filename = getunusedfilename(outputsettings.coverage.filename);
            outputsettings.mesohyperbolicity.filename = getunusedfilename(outputsettings.mesohyperbolicity.filename);
            outputsettings.gradient.filename = getunusedfilename(outputsettings.gradient.filename);
        end
    end
end

if showgrad
    gradfig = figure;
    axgrad = gca;
end

% Sample distribution to get particles and targets
[mux, muy] = sampledist(mu, nsamplepts, xlim, ylim);
[tarx, tary] = sampledist(mu, ntargets, xlim, ylim);
foundtargets = [];
findtimes = 0;
findchance = 1 - exp(-h/findtimeconst);
phi2 = [];

% Set initial positions
if size(regime, 1) > 0
    maxagents = max(regime(:, 2));
    startagents = regime(1, 2);
    n = min(numel(xti), startagents);
    xti = [xti(1:n), NaN*ones(1, maxagents - n)];
    yti = [yti(1:n), NaN*ones(1, maxagents - n)];
    geninds = isnan(xti(1:startagents)) | isnan(yti(1:startagents));
    if any(geninds)
        [xp, yp] = sampledist(mu, sum(geninds), xlim, ylim, @(n) rand(n, 1));
        xti(geninds) = xp';
        yti(geninds) = yp';
    end
else
    maxagents = numel(xti);
end
xt = [xti(1, :)*NaN; xti];
yt = [yti(1, :)*NaN; yti];

% Set up transform type and norm
[Kx, Ky] = meshgrid(0:cres-1, 0:cres-1);
Las2 = 1./(1 + Kx.^2 + Ky.^2).^(3/2);

switch lower(algorithm)
    case 'dsmc (dct, ^3/2)'
        transform = @(x) dct2(x);
        itransform = @(x) idct2(x);
        [Kx, Ky] = meshgrid(0:cres-1, 0:cres-1);
        Las = 1./(1 + Kx.^2 + Ky.^2).^(3/2);
    case 'dsmc (dct, ^1)'
        transform = @(x) dct2(x);
        itransform = @(x) idct2(x);
        [Kx, Ky] = meshgrid(0:cres-1, 0:cres-1);
        Las = 1./(1 + Kx.^2 + Ky.^2).^(2/2);
    case 'dsmc (wavelet)'
        niters = ceil(log2(cres) - 1);
        wtype = 'db16';
        bc = 'per';
        transform = @(x) fwt(fwt(x,wtype,niters,'dim',1,bc),wtype,niters,'dim',2,bc);
        itransform = @(x) ifwt(ifwt(x,wtype,niters,cres,'dim',1,bc),wtype,niters,cres,'dim',2,bc);
        [~, info] = fwt(zeros(cres),wtype,niters,'dim',1,bc);
        K = [];
        for i = 1:length(info.Lc)
            if i == 1
                scale = cres/2^(niters);
            else
                scale = cres/2^(niters+2-i);
            end
            K = [K, scale*ones(1, info.Lc(i))];
        end
        K = ones(length(K), 1)*K;
        Las = 1./(1 + K.^2 + (K').^2).^(3/2);
    case 'wavfft'
        fftscale = 4;
        niters = ceil(log2(cres) - fftscale);
        wtype = 'db16';
        bc = 'even';
        transform = @(x) wavfft2(x, wtype, niters, bc);
        itransform = @(x) iwavfft2(x, wtype, niters, bc, cres);
        [~, info] = fwt(zeros(cres),wtype,niters,'dim',1,bc);
        K = [];
        for i = 1:length(info.Lc)
            if i == 1
                scale = cres/2^(niters);
            else
                scale = cres/2^(niters+2-i);
            end
            K = [K, scale*ones(1, info.Lc(i))];
        end
        K = ones(length(K), 1)*K;
        Las = 1./(1 + K.^2 + (K').^2).^(3/2);
        K2 = (1:info.Lc(1))';
        K2 = K2 * 2^fftscale / info.Lc(1);
        K2 = K2*ones(1, info.Lc(1));
        Las(1:info.Lc(1), 1:info.Lc(1)) = 1./(1 + K2.^2 + (K2').^2).^(3/2);
    otherwise
        %error(['Unrecognized transform type: ', tftype]);
end
if length(algorithm) >= 4 && strcmp(algorithm(1:4), 'DSMC')
    algorithm = 'DSMC';
end
    

% Initialize coefficients and DSMC surface
mu = pts2img(mux, muy, cres, xlim, ylim);
logmu = max(log(max(mu, 0)/lambda), 0);
logmu = logmu / (sum(logmu(:))/cres^2);
muks = transform(logmu);

c = pts2img(xt, yt, cres, xlim, ylim);
cks = transform(c);

keepout = ones(cres);
kthick = 0;
keepout(1+kthick:end-kthick, 1+kthick:end-kthick) = 0;

sks = cks - muks;
s = itransform(Las.*sks);

% Misc setup
lastwarn('');
fcount = zeros(1, 6);
t = 0;
clear lawnmowerstep
if ~isempty(ax) || ~isempty(outputsettings)
    co = get(gca,'ColorOrder');
end

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
    [tarxn, taryn] = rk4step(t, h, tarx, tary, vx, vy);
    tdisps = sqrt((tarxn-tarx).^2 + (taryn-tary).^2);
    tarx = tarxn + tu*rand(size(tarx)).*tdisps;
    tary = taryn + tu*rand(size(tarx)).*tdisps;
    [xt, yt] = rk4step(t, h, xt, yt, vx, vy);
    
    % Get coefficients for mu
    mu = pts2img(mux, muy, cres, xlim, ylim);
    logmu = max(log(max(mu, 0)/lambda), 0);
    logmu = logmu / (sum(logmu(:))/cres^2);
    muks = transform(logmu);
    
    % Get number of agents for this step
    if isempty(regime)
        nactiveagents = numel(xti(1, :));
    else
        iregime = find(regime(:, 1) <= t, 1, 'last');
        if isempty(iregime)
            nactiveagents = 0;
        else
            nactiveagents = regime(iregime, 2);
        end
    end
    
    % Get current agent positions
    xtn0 = xt(end, 1:nactiveagents);
    ytn0 = yt(end, 1:nactiveagents);
    if any(isnan(xtn0))
        [xp, yp] = sampledist(mu, sum(isnan(xtn0)), ...
            xlim, ylim, @(n) rand(n, 1));
        xtn0(isnan(xtn0)) = xp';
        ytn0(isnan(ytn0)) = yp';
    end
    
    % Use search algorithm to determine new agent positions
    switch lower(algorithm)
        case 'dsmc'
            [xtn, ytn] = dsmcstep(xtn0, ytn0, h, umax, s, cres, xlim, ylim, au, spherical);
        case 'lawnmower'
            [xtn, ytn] = lawnmowerstep(xtn0, ytn0, h, umax, xlim, ylim, 3, spherical);
        case 'random walk'
            [xtn, ytn] = randomwalkstep(xtn0, ytn0, h, umax, xlim, ylim, spherical);
        otherwise
            error(['Unrecognized search algorithm: ', algorithm]);
    end
    npadding = maxagents - nactiveagents;
    xtn = [xtn, NaN*ones(1, npadding)];
    ytn = [ytn, NaN*ones(1, npadding)];
    xtn0 = [xtn0, NaN*ones(1, npadding)];
    ytn0 = [ytn0, NaN*ones(1, npadding)];
    xtns = bsxfun(@plus, bsxfun(@times, xtn - xtn0, linspace(0, 1, substeps + 2)'), xtn0);
    ytns = bsxfun(@plus, bsxfun(@times, ytn - ytn0, linspace(0, 1, substeps + 2)'), ytn0);
    xt = [xt; xtns(2:end, :)];
    yt = [yt; ytns(2:end, :)];
    
    % Get coefficients for coverage
    c = pts2img(xt, yt, cres, xlim, ylim);
    c2 = c + keepout*sqrt(mean(c(c > 0)));
    c2 = c2 / (sum(c2(:))/cres^2);
    cks = transform(c2);
    
    writelog(sprintf('Simulation time: %.02f / %0.2f completed.', t, maxtime));
    
    t = t + h;
    
    % Keep track of discovered targets
    if ntargets > 0
        found = findtargets(xt, yt, tarx, tary, findradius, findchance, spherical);
        found = setdiff(found, foundtargets);
        foundtargets = [foundtargets; found];
        findtimes = [findtimes; t*ones(numel(found), 1)];
    end
    
    % Calculate convergence metric
    sks = cks - muks;
    s = itransform(Las.*sks);
    sks2 = dct2(c) - dct2(logmu);
    phi2(end+1) = sum(sum(Las2.*sks2.^2));
    
    % Plot particles and trajectories
    if isfield(ax, 'main') && ~isempty(ax.main)
        plotmain(ax.main, mux, muy, tarx, tary, foundtargets, xt, yt, ...
            ntargets, findradius, xland, yland, xlim, ylim, co, spherical);
        maketitle(ax.main, t, datetitle, starttime);
    end
    
    % Plot convergence metric
    if isfield(ax, 'convergence') && ~isempty(ax.convergence)
        plotconvergence(ax.convergence, phi2);
    end
    
    % Plot log density of particles
    if isfield(ax, 'mu') && ~isempty(ax.mu)
        plotmu(ax.mu, logmu, xlim, ylim, cres);
    end
    
    % Plot trajectory density
    if isfield(ax, 'coverage') &&  ~isempty(ax.coverage)
        plotcoverage(ax.coverage, c, xlim, ylim, cres);
    end
    
    % Plot gradient
    if showgrad
        plots(axgrad, s, xlim, ylim, cres);
    end
    
    drawnow;
    
    % Save output figures
    if ~isempty(outputsettings)
        plotfun = {...
            @(ax) plotmain(ax, mux, muy, tarx, tary, foundtargets, xt, yt, ...
            ntargets, findradius, xland, yland, xlim, ylim, co, spherical), ...
            @(ax) plotconvergence(ax, phi2), ...
            @(ax) plotmu(ax, logmu, xlim, ylim, cres), ...
            @(ax) plotcoverage(ax, c, xlim, ylim, cres), ...
            @(ax) plotmesohyperbolicity(ax, vx, vy, xlim, ylim, t, ...
            outputsettings.mesohyperbolicity.T, cres)...
            @(ax) plotgrad(ax, s, xlim, ylim, cres)};
        
        osfigs = {outputsettings.main, outputsettings.convergence, ...
            outputsettings.mu, outputsettings.coverage, ...
            outputsettings.mesohyperbolicity, ...
            outputsettings.gradient};
        
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
                    [imind, map] = rgb2ind(im, 256);
                    if fcount(i) == 0
                        imwrite(imind, map, os.filename, 'DelayTime', 1/30, 'LoopCount', inf);
                    else
                        imwrite(imind, map, os.filename, 'DelayTime', 1/30, 'WriteMode', 'append');
                    end
                else
                    if ~exist([pwd, filesep, os.filename], 'file')
                        mkdir(os.filename);
                    end
                    print(outfig, [os.filename, filesep, sprintf('%.3d.png', fcount(i) + 1)], '-dpng');
                end
                fcount(i) = fcount(i) + 1;
                set(outfig, 'Visible', 'off');
            end
        end
    end
    
    % Pause
    while paused()
        pause(0.05);
    end
end

if showgrad
    close(gradfig);
end

if ~isempty(ax) || ~isempty(outputsettings)
    set(gca, 'ColorOrder', co);
end

targets = ntargets:-1:ntargets-numel(findtimes)+1;


function plotmain(ax, mux, muy, muxtar, muytar, foundtargets, xt, yt, ...
    ntargets, findradius, xland, yland, xlim, ylim, co, spherical)

plot(ax, mux, muy, '.b', 'MarkerSize', 2);
set(ax, 'ColorOrder', co(2:end, :), 'NextPlot', 'replacechildren');
hold(ax, 'on');
plot(ax, muxtar, muytar, 'co', 'MarkerSize', 5, 'LineWidth', 2);
plot(ax, muxtar(foundtargets), muytar(foundtargets), 'mo', ...
    'MarkerSize', 5, 'LineWidth', 2);

plot(ax, xland, yland, 'k.');

plot(ax, [xt(1, :); xt], [yt(1, :); yt]);
plot(ax, [xt(end, :); xt(end, :)], ...
    [yt(end, :); yt(end, :)], '.k', 'MarkerSize', 20);
plot(ax, [xt(end, :); xt(end, :)], ...
    [yt(end, :); yt(end, :)], '.', 'MarkerSize', 14);

if ntargets > 0
    theta = linspace(0, 2*pi)';
    if spherical
        circx = bsxfun(@plus, xt(end, :), ...
            bsxfun(@rdivide, findradius*cos(theta), cosd(yt(end, :))));
    else
        circx = bsxfun(@plus, xt(end, :), findradius*cos(theta));
    end
    circy = bsxfun(@plus, yt(end, :), findradius*sin(theta));
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


function plotmu(ax, mu, xlim, ylim, cres)

mu = [mu, zeros(cres, 1); zeros(1, cres), 0];
x = linspace(xlim(1), xlim(2), cres+1);
y = linspace(ylim(1), ylim(2), cres+1);
[x, y] = meshgrid(x, y);
surf(ax, x, y, mu, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Log(Particle distribution)');
caxis(ax, sort([0, max(max(mu))]));


function plotcoverage(ax, c, xlim, ylim, cres)

c = [c, zeros(cres, 1); zeros(1, cres), 0];
x = linspace(xlim(1), xlim(2), cres+1);
y = linspace(ylim(1), ylim(2), cres+1);
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


function plots(ax, s, xlim, ylim, cres)

s = [s, zeros(cres, 1); zeros(1, cres), 0];
x = linspace(xlim(1), xlim(2), cres+1);
y = linspace(ylim(1), ylim(2), cres+1);
[x, y] = meshgrid(x, y);
surf(ax, x, y, s, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'S');
title(ax, 'DSMC surface');

function plotgrad(ax, s, xlim, ylim, cres)

x = linspace(xlim(1), xlim(2), 22);
x = x(2:end-1);
y = linspace(ylim(1), ylim(2), 22);
y = y(2:end-1);
[xa, ya] = meshgrid(x, y);

Bdir = @(Bs) -bsxfun(@rdivide, Bs, sqrt(sum((Bs).^2, 1)));

tmp = Bdir(B(xa(:)', ya(:)', s, cres, xlim, ylim));
u = tmp(1, :);
v = tmp(2, :);

quiver(ax, xa(:)', ya(:)', u, v);
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Bdir');


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

if numel(filename) > 0 && exist([pwd, filesep, filename], 'file')
    i = 1;
    while exist([pwd, filesep, filename], 'file')
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

