function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% GUI for DSMC Matlab functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 1/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 09-Feb-2015 21:27:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_OpeningFcn, ...
    'gui_OutputFcn',  @gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.tab_single,'Visible','On');
set(handles.tab_params,'Visible','Off');
set(handles.tab_multiple,'Visible','Off');
set(handles.tab_output,'Visible','Off');

set(handles.panel_gaussian,'Visible','Off');
set(handles.panel_mh,'Visible','Off');

settings = struct;
data = struct;
outputsettings = struct;

settings.N = 2;
settings.xlim = [0, 1];
settings.ylim = [0, 0.5];
settings.umax = 1;
settings.tstop = 50;
settings.lambda = 0.25;
settings.algorithm = 'DSMC (DCT, ^3/2)';
settings.spherical = 0;
settings.cres = 128;
settings.h = 1/100;
settings.nsamplepts = 10000;
settings.substeps = 0;
settings.ntargets = 1000;
settings.findradius = 0.01;
settings.findtimeconst = 0.05;
settings.stopallfound = 0;
settings.mutype = 'Uniform';
settings.gaussianx = 0.5;
settings.gaussiany = 0.25;
settings.gaussianstd = 0.1;
settings.mhT = 10;
settings.mhupper = 4/(settings.mhT)^2;
settings.mhlower = 0;
settings.starttime = now;
settings.datetitle = 0;
settings.vx = @(t, x, y) -0.01*sin(2*pi*x).*cos(2*pi*y) + ...
    -0.001*cos(2*pi*t)*sin(2*pi*(x-0.25)).*cos(2*pi*(y-0.25));
settings.vy = @(t, x, y) 0.01*cos(2*pi*x).*sin(2*pi*y) + ...
    0.001*cos(2*pi*t)*cos(2*pi*(x-0.25)).*sin(2*pi*(y-0.25));
settings.agentuncertainty = 0.1;
settings.targetuncertainty = 0.1;
handles.settings = settings;
settings.mu = generatemu(handles);
[xti, yti] = sampledist(settings.mu, settings.N, settings.xlim, ...
    settings.ylim, @(n) rand(n, 1));
settings.xti = xti'; 
settings.yti = yti';
settings.regime = [];
settings.xland = [];
settings.yland = [];
settings.verbose = false;

data.findtimes = {};
data.targets = {};
data.algorithms = {};
data.convergence = {};
data.particledist = {};
data.coveragedist = {};
data.xpaths = {};
data.ypaths = {};
data.runsettings = {};

outputsettings.overwrite = 0;
outputsettings.main.enable = 0;
outputsettings.main.rate = 0;
outputsettings.main.filename = '';
outputsettings.main.animation = 0;
outputsettings.main.width = 400;
outputsettings.main.height = 400;
outputsettings.convergence = outputsettings.main;
outputsettings.mu = outputsettings.main;
outputsettings.coverage = outputsettings.main;
outputsettings.gradient = outputsettings.main;
outputsettings.mesohyperbolicity = outputsettings.main;
outputsettings.mesohyperbolicity.T = 1;

handles.settings = settings;
handles.data = data;
handles.outputsettings = outputsettings;
setparamboxes(handles);

xlabel(handles.axes_maux, 'Time');
ylabel(handles.axes_maux, 'Targets remaining');
rendermu(handles);

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setparamboxes(handles)
% Updates values in the parameter boxes in the gui

settings = handles.settings;
os = handles.outputsettings;
set(handles.edit_N, 'String', num2str(settings.N));
set(handles.edit_xlower, 'String', num2str(settings.xlim(1)));
set(handles.edit_xupper, 'String', num2str(settings.xlim(2)));
set(handles.edit_ylower, 'String', num2str(settings.ylim(1)));
set(handles.edit_yupper, 'String', num2str(settings.ylim(2)));
set(handles.edit_umax, 'String', num2str(settings.umax));
set(handles.edit_tstop, 'String', num2str(settings.tstop));
set(handles.edit_lambda, 'String', num2str(settings.lambda));
set(handles.cbox_spherical, 'Value', settings.spherical);
set(handles.edit_cres, 'String', num2str(settings.cres));
set(handles.edit_h, 'String', num2str(settings.h));
set(handles.edit_nsamplepts, 'String', num2str(settings.nsamplepts));
set(handles.edit_substeps, 'String', num2str(settings.substeps));
set(handles.edit_ntargets, 'String', num2str(settings.ntargets));
set(handles.edit_findradius, 'String', num2str(settings.findradius));
set(handles.edit_findtimeconst, 'String', num2str(settings.findtimeconst));
set(handles.edit_gaussianx, 'String', num2str(settings.gaussianx));
set(handles.edit_gaussiany, 'String', num2str(settings.gaussiany));
set(handles.edit_gaussianstd, 'String', num2str(settings.gaussianstd));
set(handles.edit_mhT, 'String', num2str(settings.mhT));
set(handles.edit_mhupper, 'String', num2str(settings.mhupper));
set(handles.edit_mhlower, 'String', num2str(settings.mhlower));
set(handles.edit_vx, 'String', char(settings.vx));
set(handles.edit_vy, 'String', char(settings.vy));
set(handles.edit_agentuncertainty, 'String', num2str(settings.agentuncertainty));
set(handles.edit_targetuncertainty, 'String', num2str(settings.targetuncertainty));
set(handles.edit_agentuncertainty, 'String', num2str(settings.agentuncertainty));
set(handles.edit_starttime, 'String', datestr(settings.starttime, 0));
set(handles.cbox_datetitle, 'Value', settings.datetitle);
set(handles.cbox_main_enable, 'Value', os.main.enable);
set(handles.cbox_main_animation, 'Value', os.main.animation);
set(handles.edit_main_rate, 'String', num2str(os.main.rate));
set(handles.edit_main_width, 'String', num2str(os.main.width));
set(handles.edit_main_height, 'String', num2str(os.main.height));
set(handles.edit_main_filename, 'String', os.main.filename);
set(handles.cbox_cov_enable, 'Value', os.coverage.enable);
set(handles.cbox_cov_animation, 'Value', os.coverage.animation);
set(handles.edit_cov_rate, 'String', num2str(os.coverage.rate));
set(handles.edit_cov_width, 'String', num2str(os.coverage.width));
set(handles.edit_cov_height, 'String', num2str(os.coverage.height));
set(handles.edit_cov_filename, 'String', os.coverage.filename);
set(handles.cbox_mu_enable, 'Value', os.mu.enable);
set(handles.cbox_mu_animation, 'Value', os.mu.animation);
set(handles.edit_mu_rate, 'String', num2str(os.mu.rate));
set(handles.edit_mu_width, 'String', num2str(os.mu.width));
set(handles.edit_mu_height, 'String', num2str(os.mu.height));
set(handles.edit_mu_filename, 'String', os.mu.filename);
set(handles.cbox_con_enable, 'Value', os.convergence.enable);
set(handles.cbox_con_animation, 'Value', os.convergence.animation);
set(handles.edit_con_rate, 'String', num2str(os.convergence.rate));
set(handles.edit_con_width, 'String', num2str(os.convergence.width));
set(handles.edit_con_height, 'String', num2str(os.convergence.height));
set(handles.edit_con_filename, 'String', os.convergence.filename);
set(handles.cbox_mh_enable, 'Value', os.mesohyperbolicity.enable);
set(handles.cbox_mh_animation, 'Value', os.mesohyperbolicity.animation);
set(handles.edit_mh_rate, 'String', num2str(os.mesohyperbolicity.rate));
set(handles.edit_mh_width, 'String', num2str(os.mesohyperbolicity.width));
set(handles.edit_mh_height, 'String', num2str(os.mesohyperbolicity.height));
set(handles.edit_mh_filename, 'String', os.mesohyperbolicity.filename);
set(handles.edit_mh_period, 'String', num2str(os.mesohyperbolicity.T));
set(handles.cbox_grad_enable, 'Value', os.gradient.enable);
set(handles.cbox_grad_animation, 'Value', os.gradient.animation);
set(handles.edit_grad_rate, 'String', num2str(os.gradient.rate));
set(handles.edit_grad_width, 'String', num2str(os.gradient.width));
set(handles.edit_grad_height, 'String', num2str(os.gradient.height));
set(handles.edit_grad_filename, 'String', os.gradient.filename);

function [vx, vy] = getvelfunction(handles)
% Get the appropriate velocity functions for expression mode and data
% interpolation mode

if get(handles.radio_expression, 'Value') == 1
    vx = handles.settings.vx;
    vy = handles.settings.vy;
elseif isfield(handles, 'veldata')
    veldata = handles.veldata;
    rearth = 6371e3; % hardcoded earth radius -- appears a few other places...
    vx = @(t, x, y) 180/pi/rearth*velocityint(t, x, y, ...
        veldata.water_u, veldata.xvals, veldata.yvals, veldata.tvals);
    vy = @(t, x, y) 180/pi/rearth./cosd(y).*velocityint(t, x, y, ...
        veldata.water_v, veldata.xvals, veldata.yvals, veldata.tvals);
else
    vx = @(x, y, t) 0;
    vy = @(x, y, t) 0;
end

function mu = generatemu(handles)
% Generates initial particle distribution according to the selected method

settings = handles.settings;
x = linspace(settings.xlim(1), settings.xlim(2), settings.cres);
y = linspace(settings.ylim(1), settings.ylim(2), settings.cres);
[x, y] = meshgrid(x, y);
if strcmp(settings.mutype, 'Uniform')
    mu = ones(size(x));
elseif strcmp(settings.mutype, 'Gaussian')
    mu = 1/settings.gaussianstd/sqrt(2*pi) * ...
        exp(-((x - settings.gaussianx).^2 + (y - settings.gaussiany).^2)/2/settings.gaussianstd^2);
elseif strcmp(settings.mutype, 'Mesohyperbolicity')
    [vx, vy] = getvelfunction(handles);
    mh = mesohyperbolicity(vx, vy, settings.xlim, ...
        settings.ylim, 0, settings.mhT, settings.cres);
    mu = zeros(size(x));
    if settings.mhlower < settings.mhupper
        mu(mh >= settings.mhlower & mh <= settings.mhupper) = 1;
    else
        mu(mh >= settings.mhlower | mh <= settings.mhupper) = 1;
    end
end
if settings.spherical
    mu = mu.*cosd(y);
end

function rendermu(handles)
% Renders a preview of the initial particle distribution

settings = handles.settings;
[mux, muy] = sampledist(settings.mu, 1000, settings.xlim, settings.ylim);
plot(handles.axes_mu, mux, muy, '.', 'MarkerSize', 2);
axis(handles.axes_mu, 'equal');
axis(handles.axes_mu, [settings.xlim, settings.ylim]);
set(handles.axes_mu, 'XTick', []);
set(handles.axes_mu, 'YTick', []);
if strcmp(settings.mutype, 'Mesohyperbolicity')
    x = linspace(settings.xlim(1), settings.xlim(2), settings.cres);
    y = linspace(settings.ylim(1), settings.ylim(2), settings.cres);
    [x, y] = meshgrid(x, y);
    [vx, vy] = getvelfunction(handles);
    mh = mesohyperbolicity(vx, vy, settings.xlim, ...
        settings.ylim, 0, settings.mhT, settings.cres);
    surf(handles.axes_mh, x, y, mh, 'EdgeColor', 'none');
    axis(handles.axes_mh, 'equal');
    axis(handles.axes_mh, [settings.xlim, settings.ylim]);
    caxis(handles.axes_mh, [-2/(settings.mhT)^2, 6/(settings.mhT)^2]);
    colorbar('peer', handles.axes_mh, 'NorthOutside');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UI button callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function btn_singlestart_Callback(hObject, ~, handles)
% Start the simulation in convergence mode

setappdata(handles.btn_singleabort, 'Abort', 0);

while ~getappdata(handles.btn_singleabort, 'Abort')
    cla(handles.axes_smain);
    cla(handles.axes_saux);
    cla(handles.axes_smu);
    cla(handles.axes_scov);
    
    settings = handles.settings;
    outputsettings = handles.outputsettings;
    data = handles.data;
    
    [vx, vy] = getvelfunction(handles);
    if settings.spherical
        rearth = 6371e3;
        umax = settings.umax/rearth*180/pi;
        findradius = settings.findradius/rearth*180/pi;
    else
        umax = settings.umax;
        findradius = settings.findradius;
    end
    
    [xti, yti] = sampledist(settings.mu, settings.N, settings.xlim, ...
        settings.ylim, @(n) rand(n, 1));
    xti = xti'; yti = yti';
    
    if isfield(handles, 'veldata')
        iland = find(isnan(handles.veldata.water_u(:,:,1)) | isnan(handles.veldata.water_v(:,:,1)));
        [xtemp, ytemp] = meshgrid(handles.veldata.xvals, handles.veldata.yvals);
        xland = xtemp(iland);
        yland = ytemp(iland);
    else
        xland = [];
        yland = [];
    end
    
    settings.ntargets = 0;
    settings.stopallfound = 0;
    settings.findradius = findradius;
    settings.umax = umax;
    settings.vx = vx;
    settings.vy = vy;
    settings.xti = xti;
    settings.yti = yti;
    settings.xland = xland;
    settings.yland = yland;
    
    ax.main = handles.axes_smain;
    ax.convergence = handles.axes_saux;
    ax.mu = handles.axes_smu;
    ax.coverage = handles.axes_scov;
    
    [~, ~, phi2, mur, c, xt, yt] = runsearch(settings, outputsettings, ...
        ax, handles.btn_singleabort, handles.cbox_spause);
    
    data.convergence = [data.convergence; {phi2}];
    data.particledist = [data.particledist; {mur}];
    data.coveragedist = [data.coveragedist; {c}];
    data.xpaths = [data.xpaths; {xt}];
    data.ypaths = [data.ypaths; {yt}];
    data.runsettings = [data.runsettings; {settings}];
    handles.data = data;
    guidata(hObject, handles);
    
    if get(handles.cbox_srunonce, 'Value')
        break;
    end
end

function btn_multistart_Callback(hObject, ~, handles)
% Start the simulation in find targets mode

setappdata(handles.btn_singleabort, 'Abort', 0);

settings = handles.settings;
outputsettings = handles.outputsettings;
data = handles.data;

[vx, vy] = getvelfunction(handles);
if settings.spherical
    rearth = 6371e3;
    umax = settings.umax/rearth*180/pi;
    findradius = settings.findradius/rearth*180/pi;
else
    umax = settings.umax;
    findradius = settings.findradius;
end

while getappdata(handles.btn_singleabort, 'Abort') == 0
    cla(handles.axes_mmain);
    
    [xti, yti] = sampledist(settings.mu, settings.N, settings.xlim, ...
        settings.ylim, @(n) rand(n, 1));
    xti = xti'; yti = yti';
    
    if isfield(handles, 'veldata')
        iland = find(isnan(handles.veldata.water_u(:,:,1)) | isnan(handles.veldata.water_v(:,:,1)));
        [xtemp, ytemp] = meshgrid(handles.veldata.xvals, handles.veldata.yvals);
        xland = xtemp(iland);
        yland = ytemp(iland);
    else
        xland = [];
        yland = [];
    end
    
    settings.stopallfound = 1;
    settings.findradius = findradius;
    settings.umax = umax;
    settings.vx = vx;
    settings.vy = vy;
    settings.xti = xti;
    settings.yti = yti;
    settings.xland = xland;
    settings.yland = yland;
    
    ax.main = handles.axes_mmain;
    
    [time, target] = runsearch(settings, outputsettings, ...
        ax, handles.btn_singleabort, handles.cbox_mpause);
    
    hold(handles.axes_maux, 'on');
    if strcmp(settings.algorithm, 'DSMC')
        plot(handles.axes_maux, time, target, 'b.');
    elseif strcmp(settings.algorithm, 'Lawnmower')
        plot(handles.axes_maux, time, target, 'rx');
    else % random walk
        plot(handles.axes_maux, time, target, 'go');
    end
    hold(handles.axes_maux, 'off');
    
    data.findtimes = [data.findtimes; {time}];
    data.targets = [data.targets; {target}];
    data.algorithms = [data.algorithms; {settings.algorithm}];
    handles.data = data;
    guidata(hObject, handles);
end

function btn_save_Callback(~, ~, handles)

settings = handles.settings;
data = handles.data;
uisave({'data', 'settings'});

function btn_load_Callback(hObject, ~, handles)

uiopen('load');
if exist('data', 'var') && exist('settings', 'var')
    if isfield(settings, 'L1') && isfield(settings, 'L2')
        settings.xlim = [0, settings.L1];
        settings.ylim = [0, settings.L2];
        settings = rmfield(settings, {'L1', 'L2'});
    end 
    if ~isfield(settings, 'starttime')
        settings.starttime = now;
        settings.datetitle = 0;
    end
    if isfield(settings, 'gifname')
        settings = rmfield(settings, 'gifname');
    end
    
    handles.settings = settings;
    handles.data = data;
    
    if strcmp(settings.mutype, 'Gaussian')
        set(handles.panel_gaussian,'Visible','On');
        set(handles.panel_mh,'Visible','Off');
    elseif strcmp(settings.mutype, 'Mesohyperbolicity');
        set(handles.panel_gaussian,'Visible','Off');
        set(handles.panel_mh,'Visible','On');
    else
        set(handles.panel_gaussian,'Visible','Off');
        set(handles.panel_mh,'Visible','Off');
    end
    
    rendermu(handles);
    
    cla(handles.axes_maux);
    n = numel(data.findtimes);
    hold(handles.axes_maux, 'on');
    for i = 1:n
        if strcmp(data.algorithms{i}, 'DSMC')
            plot(handles.axes_maux, data.findtimes{i}, data.targets{i}, 'b.');
        elseif strcmp(data.algorithms{i}, 'Lawnmower')
            plot(handles.axes_maux, data.findtimes{i}, data.targets{i}, 'rx');
        else % random walk
            plot(handles.axes_maux, data.findtimes{i}, data.targets{i}, 'go');
        end
    end
    hold(handles.axes_maux, 'off');
    xlabel(handles.axes_maux, 'Time');
    ylabel(handles.axes_maux, 'Targets remaining');
    
    setparamboxes(handles);
    
    guidata(hObject, handles);
end

function btn_mclear_Callback(hObject, ~, handles)

button = questdlg('Clear data from previous runs; are you sure?','Confirm Clear','Cancel');

if strcmp(button, 'Yes')
    cla(handles.axes_maux);
    handles.data.findtimes = {};
    handles.data.targets = {};
    handles.data.algorithms = {};
    xlabel(handles.axes_maux, 'Time');
    ylabel(handles.axes_maux, 'Targets remaining');
    guidata(hObject, handles);
end

function btn_loadveldata_Callback(hObject, ~, handles)

[filename, path] = uigetfile('*.mat');
if filename
    h = waitbar(0, 'Loading velocity data...');
    load([path, filename]);
    if exist('veldata', 'var')
        waitbar(1, h, 'Loaded velocity data.');
        handles.veldata = veldata;
        set(handles.text_velfile, 'String', filename);
        guidata(hObject, handles);
        pause(0.5);
        close(h);
    else
        waitbar(1, h, 'Failed to load velocity data.');
        pause(0.5);
        close(h);
    end
end

function uipanel3_SelectionChangeFcn(~, eventdata, handles)
% Handle tab changes and change the colormap -- only one colormap can be
% associated with a single figure at once, despite having multiple axes

switch eventdata.NewValue
    case handles.tabradio_params
        set(handles.tab_params,'Visible','On');
        set(handles.tab_single,'Visible','Off');
        set(handles.tab_multiple,'Visible','Off');
        set(handles.tab_output,'Visible','Off');
        colormap(mhcolormap(256));
    case handles.tabradio_single
        set(handles.tab_params,'Visible','Off');
        set(handles.tab_single,'Visible','On');
        set(handles.tab_multiple,'Visible','Off');
        set(handles.tab_output,'Visible','Off');
        colormap(jet(256));
    case handles.tabradio_multiple
        set(handles.tab_params,'Visible','Off');
        set(handles.tab_single,'Visible','Off');
        set(handles.tab_multiple,'Visible','On');
        set(handles.tab_output,'Visible','Off');
        colormap(jet(256));
    case handles.tabradio_output
        set(handles.tab_params,'Visible','Off');
        set(handles.tab_single,'Visible','Off');
        set(handles.tab_multiple,'Visible','Off');
        set(handles.tab_output,'Visible','On');
        colormap(jet(256));
end

function uipanel_velradio_SelectionChangeFcn(~, eventdata, handles)

switch eventdata.NewValue
    case handles.radio_expression
        set(handles.edit_vx,'Enable','On');
        set(handles.edit_vy,'Enable','On');
    case handles.radio_interpolation
        set(handles.edit_vx,'Enable','Off');
        set(handles.edit_vy,'Enable','Off');
end

function btn_singleabort_Callback(~, ~, handles)

setappdata(handles.btn_singleabort, 'Abort', 1);

function btn_multiabort_Callback(~, ~, handles)

setappdata(handles.btn_singleabort, 'Abort', 1); % not a bug

function btn_sresume_Callback(~, ~, ~)
% Blank callback function

function cbox_spause_Callback(~, ~, ~)
% Blank callback function

function cbox_mpause_Callback(~, ~, ~)
% Blank callback function

function cbox_srunonce_Callback(~, ~, ~)
% Blank callback function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter edit callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_N_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.N = val;
guidata(hObject, handles);

function edit_umax_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.umax = val;
guidata(hObject, handles);

function edit_tstop_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.tstop = val;
guidata(hObject, handles);

function edit_ntargets_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.ntargets = val;
guidata(hObject, handles);

function edit_findradius_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.findradius = val;
guidata(hObject, handles);

function edit_findtimeconst_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.findtimeconst = val;
guidata(hObject, handles);

function edit_cres_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.cres = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_h_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.h = val;
guidata(hObject, handles);

function edit_nsamplepts_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.nsamplepts = val;
guidata(hObject, handles);

function edit_substeps_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = floor(eval(str));
set(hObject, 'String', num2str(val));
handles.settings.substeps = val;
guidata(hObject, handles);

function edit_lambda_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.lambda = val;
guidata(hObject, handles);

function popup_algorithm_Callback(hObject, ~, handles)

contents = cellstr(get(hObject,'String'));
str = contents{get(hObject,'Value')};
handles.settings.algorithm = str;
guidata(hObject, handles);

function edit_vx_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', char(val));
handles.settings.vx = val;
if strcmp(handles.settings.mutype, 'Mesohyperbolicity')
    handles.settings.mu = generatemu(handles);
    rendermu(handles);
end
guidata(hObject, handles);

function edit_vy_Callback(hObject, ~, handles) %#ok<*DEFNU>

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', char(val));
handles.settings.vy = val;
if strcmp(handles.settings.mutype, 'Mesohyperbolicity')
    handles.settings.mu = generatemu(handles);
    rendermu(handles);
end
guidata(hObject, handles);

function edit_agentuncertainty_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.agentuncertainty = val;
guidata(hObject, handles);

function edit_targetuncertainty_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.targetuncertainty = val;
guidata(hObject, handles);

function popup_mutype_Callback(hObject, ~, handles)

contents = cellstr(get(hObject,'String'));
str = contents{get(hObject,'Value')};
handles.settings.mutype = str;
handles.settings.mu = generatemu(handles);
rendermu(handles);
settings = handles.settings;
if strcmp(settings.mutype, 'Gaussian')
    set(handles.panel_gaussian,'Visible','On');
    set(handles.panel_mh,'Visible','Off');
elseif strcmp(settings.mutype, 'Mesohyperbolicity');
    set(handles.panel_gaussian,'Visible','Off');
    set(handles.panel_mh,'Visible','On');
else
    set(handles.panel_gaussian,'Visible','Off');
    set(handles.panel_mh,'Visible','Off');
end
guidata(hObject, handles);

function edit_gaussianstd_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussianstd = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_gaussianx_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussianx = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_gaussiany_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussiany = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_mhlower_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhlower = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_mhupper_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhupper = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_mhT_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhT = val;
handles.settings.mu = generatemu(handles);
rendermu(handles);
guidata(hObject, handles);

function edit_xupper_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
%if val > handles.settings.xlim(1)
    set(hObject, 'String', num2str(val));
    handles.settings.xlim(2) = val;
    handles.settings.mu = generatemu(handles);
    guidata(hObject, handles);
    rendermu(handles);
%end

function edit_xlower_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
%if val < handles.settings.xlim(2)
    set(hObject, 'String', num2str(val));
    handles.settings.xlim(1) = val;
    handles.settings.mu = generatemu(handles);
    guidata(hObject, handles);
    rendermu(handles);
%end

function edit_yupper_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
%if val > handles.settings.ylim(1)
    set(hObject, 'String', num2str(val));
    handles.settings.ylim(2) = val;
    handles.settings.mu = generatemu(handles);
    guidata(hObject, handles);
    rendermu(handles);
%end

function edit_ylower_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
%if val < handles.settings.ylim(2)
    set(hObject, 'String', num2str(val));
    handles.settings.ylim(1) = val;
    handles.settings.mu = generatemu(handles);
    guidata(hObject, handles);
    rendermu(handles);
%end

function cbox_spherical_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.settings.spherical = val;
guidata(hObject, handles);

function cbox_datetitle_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.settings.datetitle = val;
guidata(hObject, handles);

function edit_starttime_Callback(hObject, ~, handles)

val = datenum(get(hObject, 'String'));
set(hObject, 'String', datestr(val, 0));
handles.settings.starttime = val;
guidata(hObject, handles);

function cbox_overwrite_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.overwrite = val;
guidata(hObject, handles);

function cbox_main_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.main.enable = val;
guidata(hObject, handles);

function edit_main_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.rate = val;
guidata(hObject, handles);

function edit_main_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.width = val;
guidata(hObject, handles);

function edit_main_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.height = val;
guidata(hObject, handles);

function edit_main_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.main.filename = str;
guidata(hObject, handles);

function cbox_main_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.main.animation = val;
guidata(hObject, handles);

function cbox_cov_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.coverage.enable = val;
guidata(hObject, handles);

function edit_cov_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.rate = val;
guidata(hObject, handles);

function edit_cov_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.width = val;
guidata(hObject, handles);

function edit_cov_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.height = val;
guidata(hObject, handles);

function edit_cov_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.coverage.filename = str;
guidata(hObject, handles);

function cbox_cov_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.coverage.animation = val;
guidata(hObject, handles);

function cbox_mu_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.mu.enable = val;
guidata(hObject, handles);

function edit_mu_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.rate = val;
guidata(hObject, handles);

function edit_mu_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.width = val;
guidata(hObject, handles);

function edit_mu_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.height = val;
guidata(hObject, handles);

function edit_mu_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.mu.filename = str;
guidata(hObject, handles);

function cbox_mu_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.mu.animation = val;
guidata(hObject, handles);

function cbox_con_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.convergence.enable = val;
guidata(hObject, handles);

function edit_con_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.rate = val;
guidata(hObject, handles);

function edit_con_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.width = val;
guidata(hObject, handles);

function edit_con_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.height = val;
guidata(hObject, handles);

function edit_con_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.convergence.filename = str;
guidata(hObject, handles);

function cbox_con_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.convergence.animation = val;
guidata(hObject, handles);

function cbox_mh_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.mesohyperbolicity.enable = val;
guidata(hObject, handles);

function edit_mh_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.rate = val;
guidata(hObject, handles);

function edit_mh_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.width = val;
guidata(hObject, handles);

function edit_mh_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.height = val;
guidata(hObject, handles);

function edit_mh_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.mesohyperbolicity.filename = str;
guidata(hObject, handles);

function cbox_mh_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.mesohyperbolicity.animation = val;
guidata(hObject, handles);

function edit_mh_period_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.T = val;
guidata(hObject, handles);

function cbox_grad_enable_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.gradient.enable = val;
guidata(hObject, handles);

function edit_grad_rate_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.grady.rate = val;
guidata(hObject, handles);

function edit_grad_width_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.gradient.width = val;
guidata(hObject, handles);

function edit_grad_height_Callback(hObject, ~, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.gradient.height = val;
guidata(hObject, handles);

function edit_grad_filename_Callback(hObject, ~, handles)

str = get(hObject,'String');
handles.outputsettings.gradient.filename = str;
guidata(hObject, handles);

function cbox_grad_animation_Callback(hObject, ~, handles)

val = get(hObject, 'Value');
handles.outputsettings.gradient.animation = val;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UI element creation functions, all boilerplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_N_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_umax_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_tstop_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ntargets_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_findradius_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_findtimeconst_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cres_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_h_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_nsamplepts_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_lambda_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popup_algorithm_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_vx_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_vy_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_agentuncertainty_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_targetuncertainty_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popup_mutype_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussianstd_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussianx_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussiany_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhlower_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhupper_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhT_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_xupper_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_xlower_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_yupper_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ylower_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_starttime_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_rate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_width_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_height_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_filename_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_rate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_width_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_height_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_filename_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_rate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_width_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_height_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_filename_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_rate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_width_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_height_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_filename_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_filename_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_height_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_width_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_rate_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_period_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_substeps_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_grad_rate_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_grad_width_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_grad_height_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_grad_filename_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
