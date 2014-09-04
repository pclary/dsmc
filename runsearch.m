function [findtimes, targets, phi2, mur, c, xt1, xt2] = runsearch(settings, ...
    outputsettings, ax, aborthandle, pausebutton)
%RUNSEARCH Runs a search using the specified algorithm
%   Returns the time at which each target was found and a vector of the
%   number of remaining targets at each detection time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 7/23/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack settings
xlim = settings.xlim;
ylim = settings.ylim;
K = settings.K;
v1 = settings.v1;
v2 = settings.v2;
h = settings.h;
x1 = settings.x1;
x2 = settings.x2;
mu = settings.mu;
nsamplepts = settings.nsamplepts;
ngrid = settings.ngrid;
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
xt1i = settings.xt1i;
xt2i = settings.xt2i;
xland = settings.xland;
yland = settings.yland;

% Append (1), (2), etc to output file names to avoid overwriting old output
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

% Precompute Fourier basis functions
fks = cell(K + 1);
hks = zeros(K + 1);
for K1 = 0:K
    for K2 = 0:K
        k1 = K1*pi/(xlim(2)-xlim(1));
        k2 = K2*pi/(ylim(2)-ylim(1));
        hks(K1+1, K2+1) = sqrt(trapz(x2(:,1), trapz(x1(1,:), cos(k1*(x1-xlim(1))).^2.*cos(k2*(x2-ylim(1))).^2, 2), 1));
        fks{K1+1, K2+1} = 1./hks(K1+1, K2+1).*cos(k1*(x1-xlim(1))).*cos(k2*(x2-ylim(1)));
    end
end

% Sample distribution to get particles and targets
[mu1, mu2] = sampledist(mu, x1(1, :), x2(:, 1), nsamplepts);
[mu1tar, mu2tar] = sampledist(mu, x1(1, :), x2(:, 1), ntargets);
foundtargets = [];
findtimes = 0;
findchance = 1 - exp(-h/findtimeconst);
phi2 = [];

% Set initial positions
xt1 = xt1i;
xt2 = xt2i;
if strcmp(algorithm, 'Lawnmower')
    xt1c = xt1i;
    xt2c = xt2i;
end

fcount = [0, 0, 0, 0, 0];
stepnum = size(xt1, 1);
t = 0;
co = get(gca,'ColorOrder');

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
    [mu1, mu2] = rk4step(t, h, mu1, mu2, v1, v2);
    [mu1tarn, mu2tarn] = rk4step(t, h, mu1tar, mu2tar, v1, v2);
    tdisps = sqrt((mu1tarn-mu1tar).^2 + (mu2tarn-mu2tar).^2);
    mu1tar = mu1tarn + tu*rand(size(mu1tar)).*tdisps;
    mu2tar = mu2tarn + tu*rand(size(mu1tar)).*tdisps;
    [xt1, xt2] = rk4step(t, h, xt1, xt2, v1, v2);
    
    % Use search algorithm to determine new agent positions
    if strcmp(algorithm, 'DSMC')
        xt10 = xt1(1:stepnum, :);
        xt20 = xt2(1:stepnum, :);
        muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        [xt1n, xt2n] = dsmcstep(xt10, xt20, muks, h, umax, hks, K, xlim, ylim, au, spherical);
    elseif strcmp(algorithm, 'Lawnmower')
        % Keep track of unadvected trajectory (required by lawnmowerstep)
        xt10 = xt1c(1:stepnum, :);
        xt20 = xt2c(1:stepnum, :);
        xt10(end, :) = xt1(end, :);
        xt20(end, :) = xt2(end, :);
        [xt1n, xt2n] = lawnmowerstep(xt10, xt20, h, umax, xlim, ylim, 10, spherical);
        xt1c(stepnum+1, :) = xt1n;
        xt2c(stepnum+1, :) = xt2n;
        % Update mu only every 10th step (significant speedup)
        if mod(stepnum, 10) == 1
            muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        end
    elseif strcmp(algorithm, 'Random Walk')
        xt10 = xt1(1:stepnum, :);
        xt20 = xt2(1:stepnum, :);
        [xt1n, xt2n] = randomwalkstep(xt10, xt20, h, umax, xlim, ylim, spherical);
        % Update mu only every 10th step (significant speedup)
        if mod(stepnum, 10) == 1
            muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks);
        end
    end
    xt1(stepnum+1, :) = xt1n;
    xt2(stepnum+1, :) = xt2n;
    stepnum = stepnum + 1;
    
    t = t + h;
    
    % Keep track of discovered targets
    notfound = ~ismember(1:ntargets, foundtargets);
    dists = min(sqrt(bsxfun(@minus, xt1(stepnum, :), mu1tar).^2 + ...
        bsxfun(@minus, xt2(stepnum, :), mu2tar).^2), [], 2);
    if numel(dists) == 0
        dists = zeros(0, 1);
    end
    found = find(dists < findradius & notfound' & rand(ntargets, 1) < findchance);
    foundtargets = [foundtargets; found];
    findtimes = [findtimes; t*ones(numel(found), 1)];
    
    % Calculate convergence metric
    sk = @(K1, K2) ck(K1, K2, xt1(1:stepnum, :), xt2(1:stepnum, :), ...
        hks(K1+1, K2+1), xlim, ylim) - muks(K1+1, K2+1);
    phi = 0;
    for K1 = 0:K
        for K2 = 0:K
            La = 1/(1 + K1^2 + K2^2)^(3/2);
            phi = phi + La*sk(K1, K2)^2;
        end
    end
    phi2 = [phi2; phi];
    
    
    % Plot particles and trajectories
    if isfield(ax, 'main') && ~isempty(ax.main)
        plotmain(ax.main, mu1, mu2, mu1tar, mu2tar, foundtargets, xt1, xt2, ...
            stepnum, ntargets, findradius, xland, yland, xlim, ylim, co);
        maketitle(ax.main, t, datetitle, starttime);
    end
    
    % Plot convergence metric
    if isfield(ax, 'convergence') && ~isempty(ax.convergence)
        plotconvergence(ax.convergence, phi2);
    end
    
    % Plot log density of particles
    if isfield(ax, 'mu') && ~isempty(ax.mu)
        plotmu(ax.mu, mu, fks, muks, K, x1, x2, xlim, ylim);
    end
    
    % Plot trajectory density
    if isfield(ax, 'coverage') &&  ~isempty(ax.coverage)
        plotcoverage(ax.coverage, mu, fks, hks, K, xt1, xt2, x1, x2, xlim, ylim, stepnum);
    end
    
    drawnow;
    
    % Save output figures
    plotfun = {...
        @(ax) plotmain(ax, mu1, mu2, mu1tar, mu2tar, foundtargets, xt1, xt2, ...
        stepnum, ntargets, findradius, xland, yland, xlim, ylim, co), ...
        @(ax) plotconvergence(ax, phi2), ...
        @(ax) plotmu(ax, mu, fks, muks, K, x1, x2, xlim, ylim), ...
        @(ax) plotcoverage(ax, mu, fks, hks, K, xt1, xt2, x1, x2, xlim, ylim, stepnum), ...
        @(ax) plotmesohyperbolicity(ax, v1, v2, xlim, ylim, t, outputsettings.mesohyperbolicity.T, ngrid)};
    
    osfigs = {outputsettings.main, outputsettings.convergence, ...
        outputsettings.mu, outputsettings.coverage, ...
        outputsettings.mesohyperbolicity};
    
    f = figure;
    axnew = gca;
    set (f, 'Renderer', 'zbuffer');
    set(f, 'Visible', 'off');
    
    for i = 1:numel(osfigs)
        os = osfigs{i};
        if os.enable && t >= fcount(i)*os.rate
            set(f, 'Position', [0, 0, os.width, os.height]);
            func = plotfun{i};
            func(axnew);
            maketitle(axnew, t, datetitle, starttime);
            set(f, 'Visible', 'off');
            drawnow;
            if os.animation
                frame = getframe(f);
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
                frame = getframe(f);
                im = frame2im(frame);
                imwrite(im, [os.filename, '/', sprintf('%.3d.png', fcount(i) + 1)], 'png');
            end
            fcount(i) = fcount(i) + 1;
        end
    end
    delete(f);
    
    % Pause
    while paused()
        pause(0.05);
    end
end

set(gca, 'ColorOrder', co);

targets = ntargets:-1:ntargets-numel(findtimes)+1;

mur = zeros(size(mu));
for K1 = 0:K
    for K2 = 0:K
        mur = mur + fks{K1+1, K2+1}*muks(K1+1, K2+1);
    end
end
mur = mur / trapz(x2(:, 1), trapz(x1(1, :), mur, 2));

c = zeros(size(mu));
for K1 = 0:K
    for K2 = 0:K
        c = c + fks{K1+1, K2+1}*ck(K1, K2, xt1(1:stepnum, :), ...
            xt2(1:stepnum, :), hks(K1+1, K2+1), xlim, ylim);
    end
end
c = c / trapz(x2(:, 1), trapz(x1(1, :), c, 2));


function plotmain(ax, mu1, mu2, mu1tar, mu2tar, foundtargets, xt1, xt2, ...
    stepnum, ntargets, findradius, xland, yland, xlim, ylim, co)

plot(ax, mu1, mu2, '.b', 'MarkerSize', 2);
set(ax, 'ColorOrder', co(2:end, :), 'NextPlot', 'replacechildren');
hold(ax, 'on');
plot(ax, mu1tar, mu2tar, 'co', 'MarkerSize', 5, 'LineWidth', 2);
plot(ax, mu1tar(foundtargets), mu2tar(foundtargets), 'mo', ...
    'MarkerSize', 5, 'LineWidth', 2);

plot(ax, xland, yland, 'k.');

plot(ax, [xt1(1, :); xt1(1:stepnum, :)], [xt2(1, :); xt2(1:stepnum, :)]);
plot(ax, [xt1(stepnum, :); xt1(stepnum, :)], ...
    [xt2(stepnum, :); xt2(stepnum, :)], '.', 'MarkerSize', 15);

if ntargets > 0
    theta = linspace(0, 2*pi)';
    circx = bsxfun(@plus, xt1(stepnum, :), findradius*cos(theta));
    circy = bsxfun(@plus, xt2(stepnum, :), findradius*sin(theta));
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


function plotmu(ax, mu, fks, muks, K, x1, x2, xlim, ylim)

mur = zeros(size(mu));
for K1 = 0:K
    for K2 = 0:K
        mur = mur + fks{K1+1, K2+1}*muks(K1+1, K2+1);
    end
end
mur = mur / trapz(x2(:, 1), trapz(x1(1, :), mur, 2));
surf(ax, x1, x2, mur, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Log(Particle distribution)');
caxis(ax, sort([0, max(max(mur))]));


function plotcoverage(ax, mu, fks, hks, K, xt1, xt2, x1, x2, xlim, ylim, stepnum)

c = zeros(size(mu));
for K1 = 0:K
    for K2 = 0:K
        c = c + fks{K1+1, K2+1}*ck(K1, K2, xt1(1:stepnum, :), ...
            xt2(1:stepnum, :), hks(K1+1, K2+1), xlim, ylim);
    end
end
c = c / trapz(x2(:, 1), trapz(x1(1, :), c, 2));
surf(ax, x1, x2, c, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Coverage density');
maxc = max(max(c));
if isnan(maxc)
    maxc = 1;
end
caxis(ax, [0, maxc]);


function plotmesohyperbolicity(ax, vx, vy, xlim, ylim, t, T, ngrid)

[mh, x0, y0] = mesohyperbolicity(vx, vy, xlim, ylim, t, T, ngrid);
surf(ax, x0, y0, mh, 'EdgeColor', 'none');
axis(ax, 'equal');
axis(ax, [xlim, ylim]);
title(ax, 'Mesohyperbolicity');
colormap(ax, mhcolormap(256));
caxis(ax, [-2/(T)^2, 6/(T)^2]);
colorbar('peer', ax, 'EastOutside');


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

