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
% Updated 7/11/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 11-Jul-2014 00:47:26

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
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
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
settings.lambda = 0.5;
settings.algorithm = 'DSMC';
settings.spherical = 0;
settings.K = 20;
settings.ngrid = 40;
settings.h = 1/30;
settings.nsamplepts = 10000;
settings.ntargets = 100;
settings.findradius = 0.03;
settings.findtimeconst = 0.02;
settings.mutype = 'Uniform';
settings.gaussianx = 0.5;
settings.gaussiany = 0.25;
settings.gaussianstd = 0.1;
settings.mhT = 10;
settings.mhupper = 4/(settings.mhT)^2;
settings.mhlower = 0;
settings.starttime = now;
settings.datetitle = 0;
settings.v1 = @(t, x, y) -0.1*sin(2*pi*x).*cos(2*pi*y) + ...
    -0.01*cos(2*pi*t)*sin(2*pi*(x-0.25)).*cos(2*pi*(y-0.25));
settings.v2 = @(t, x, y) 0.1*cos(2*pi*x).*sin(2*pi*y) + ...
    0.01*cos(2*pi*t)*cos(2*pi*(x-0.25)).*sin(2*pi*(y-0.25));
settings.agentuncertainty = 0.1;
settings.targetuncertainty = 0.1;
handles.settings = settings;
[settings.mu, settings.x1, settings.x2] = generatemu(handles);

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
function varargout = gui_OutputFcn(hObject, eventdata, handles)
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
set(handles.edit_K, 'String', num2str(settings.K));
set(handles.edit_ngrid, 'String', num2str(settings.ngrid));
set(handles.edit_h, 'String', num2str(settings.h));
set(handles.edit_nsamplepts, 'String', num2str(settings.nsamplepts));
set(handles.edit_ntargets, 'String', num2str(settings.ntargets));
set(handles.edit_findradius, 'String', num2str(settings.findradius));
set(handles.edit_findtimeconst, 'String', num2str(settings.findtimeconst));
set(handles.edit_gaussianx, 'String', num2str(settings.gaussianx));
set(handles.edit_gaussiany, 'String', num2str(settings.gaussiany));
set(handles.edit_gaussianstd, 'String', num2str(settings.gaussianstd));
set(handles.edit_mhT, 'String', num2str(settings.mhT));
set(handles.edit_mhupper, 'String', num2str(settings.mhupper));
set(handles.edit_mhlower, 'String', num2str(settings.mhlower));
set(handles.edit_v1, 'String', char(settings.v1));
set(handles.edit_v2, 'String', char(settings.v2));
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

function [vx, vy] = getvelfunction(handles)
% Get the appropriate velocity functions for expression mode and data
% interpolation mode

if get(handles.radio_expression, 'Value') == 1
    vx = handles.settings.v1;
    vy = handles.settings.v2;
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

function [mu, x1, x2] = generatemu(handles)
% Generates initial particle distribution according to the selected method

settings = handles.settings;
x1 = linspace(settings.xlim(1), settings.xlim(2), settings.ngrid);
x2 = linspace(settings.ylim(1), settings.ylim(2), settings.ngrid);
[x1, x2] = meshgrid(x1, x2);
if strcmp(settings.mutype, 'Uniform')
    mu = ones(size(x1));
elseif strcmp(settings.mutype, 'Gaussian')
    mu = 1/settings.gaussianstd/sqrt(2*pi) * ...
        exp(-((x1 - settings.gaussianx).^2 + (x2 - settings.gaussiany).^2)/2/settings.gaussianstd^2);
elseif strcmp(settings.mutype, 'Mesohyperbolicity')
    [vx, vy] = getvelfunction(handles);
    mh = mesohyperbolicity(vx, vy, settings.xlim, ...
        settings.ylim, 0, settings.mhT, settings.ngrid);
    mu = zeros(size(x1));
    if settings.mhlower < settings.mhupper
        mu(mh >= settings.mhlower & mh <= settings.mhupper) = 1;
    else
        mu(mh >= settings.mhlower | mh <= settings.mhupper) = 1;
    end
end

function rendermu(handles)
% Renders a preview of the initial particle distribution

settings = handles.settings;
x1 = linspace(settings.xlim(1), settings.xlim(2), settings.ngrid);
x2 = linspace(settings.ylim(1), settings.ylim(2), settings.ngrid);
[mu1, mu2] = sampledist(settings.mu, x1, x2, 1000);
plot(handles.axes_mu, mu1, mu2, '.', 'MarkerSize', 2);
axis(handles.axes_mu, 'equal');
axis(handles.axes_mu, [settings.xlim, settings.ylim]);
set(handles.axes_mu, 'XTick', []);
set(handles.axes_mu, 'YTick', []);
if strcmp(settings.mutype, 'Mesohyperbolicity')
    [x1, x2] = meshgrid(x1, x2);
    [vx, vy] = getvelfunction(handles);
    mh = mesohyperbolicity(vx, vy, settings.xlim, ...
        settings.ylim, 0, settings.mhT, settings.ngrid);
    surf(handles.axes_mh, x1, x2, mh, 'EdgeColor', 'none');
    axis(handles.axes_mh, 'equal');
    axis(handles.axes_mh, [settings.xlim, settings.ylim]);
    caxis(handles.axes_mh, [-2/(settings.mhT)^2, 6/(settings.mhT)^2]);
    colorbar('peer', handles.axes_mh, 'NorthOutside');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UI button callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function btn_singlestart_Callback(hObject, eventdata, handles)
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
    
    x1 = linspace(settings.xlim(1), settings.xlim(2), settings.ngrid);
    x2 = linspace(settings.ylim(1), settings.ylim(2), settings.ngrid);
    [xt1i, xt2i] = sampledist(settings.mu, x1, x2, settings.N, @(n) rand(n, 1));
    xt1i = xt1i'; xt2i = xt2i';
    
    settings.ntargets = 0;
    settings.stopallfound = 0;
    settings.findradius = findradius;
    settings.umax = umax;
    settings.v1 = vx;
    settings.v2 = vy;
    settings.xt1i = xt1i;
    settings.xt2i = xt2i;
    
    ax.main = handles.axes_smain;
    ax.convergence = handles.axes_saux;
    ax.mu = handles.axes_smu;
    ax.coverage = handles.axes_scov;
    
    [~, ~, phi2, mur, c, xt1, xt2] = runsearch(settings, outputsettings, ...
        ax, handles.btn_singleabort, handles.cbox_spause);
    
    data.convergence = [data.convergence; {phi2}];
    data.particledist = [data.particledist; {mur}];
    data.coveragedist = [data.coveragedist; {c}];
    data.xpaths = [data.xpaths; {xt1}];
    data.ypaths = [data.ypaths; {xt2}];
    data.runsettings = [data.runsettings; {settings}];
    handles.data = data;
    guidata(hObject, handles);
    
    if get(handles.cbox_srunonce, 'Value')
        break;
    end
end

function btn_multistart_Callback(hObject, eventdata, handles)
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
    
    x1 = linspace(settings.xlim(1), settings.xlim(2), settings.ngrid);
    x2 = linspace(settings.ylim(1), settings.ylim(2), settings.ngrid);
    [xt1i, xt2i] = sampledist(settings.mu, x1, x2, settings.N, @(n) rand(n, 1));
    xt1i = xt1i'; xt2i = xt2i';
    
    settings.stopallfound = 1;
    settings.findradius = findradius;
    settings.umax = umax;
    settings.v1 = vx;
    settings.v2 = vy;
    settings.xt1i = xt1i;
    settings.xt2i = xt2i;
    
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

function btn_save_Callback(hObject, eventdata, handles)

settings = handles.settings;
data = handles.data;
uisave({'data', 'settings'});

function btn_load_Callback(hObject, eventdata, handles)

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

function btn_mclear_Callback(hObject, eventdata, handles)

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

function btn_loadveldata_Callback(hObject, eventdata, handles)

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

function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
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

function uipanel_velradio_SelectionChangeFcn(hObject, eventdata, handles)

switch eventdata.NewValue
    case handles.radio_expression
        set(handles.edit_v1,'Enable','On');
        set(handles.edit_v2,'Enable','On');
    case handles.radio_interpolation
        set(handles.edit_v1,'Enable','Off');
        set(handles.edit_v2,'Enable','Off');
end

function btn_singleabort_Callback(hObject, eventdata, handles)

setappdata(handles.btn_singleabort, 'Abort', 1);

function btn_multiabort_Callback(hObject, eventdata, handles)

setappdata(handles.btn_singleabort, 'Abort', 1); % not a bug

function btn_sresume_Callback(hObject, eventdata, handles)
% Blank callback function

function cbox_spause_Callback(hObject, eventdata, handles)
% Blank callback function

function cbox_mpause_Callback(hObject, eventdata, handles)
% Blank callback function

function cbox_srunonce_Callback(hObject, eventdata, handles)
% Blank callback function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter edit callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_N_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.N = val;
guidata(hObject, handles);

function edit_umax_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.umax = val;
guidata(hObject, handles);

function edit_tstop_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.tstop = val;
guidata(hObject, handles);

function edit_ntargets_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.ntargets = val;
guidata(hObject, handles);

function edit_K_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.K = val;
guidata(hObject, handles);

function edit_findradius_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.findradius = val;
guidata(hObject, handles);

function edit_findtimeconst_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.findtimeconst = val;
guidata(hObject, handles);

function edit_ngrid_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.ngrid = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_h_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.h = val;
guidata(hObject, handles);

function edit_nsamplepts_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.nsamplepts = val;
guidata(hObject, handles);

function edit_lambda_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.lambda = val;
guidata(hObject, handles);

function popup_algorithm_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'));
str = contents{get(hObject,'Value')};
handles.settings.algorithm = str;
guidata(hObject, handles);

function edit_v1_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', char(val));
handles.settings.v1 = val;
if strcmp(handles.settings.mutype, 'Mesohyperbolicity')
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
end
guidata(hObject, handles);

function edit_v2_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', char(val));
handles.settings.v2 = val;
if strcmp(handles.settings.mutype, 'Mesohyperbolicity')
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
end
guidata(hObject, handles);

function edit_agentuncertainty_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.agentuncertainty = val;
guidata(hObject, handles);

function edit_targetuncertainty_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.targetuncertainty = val;
guidata(hObject, handles);

function popup_mutype_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'));
str = contents{get(hObject,'Value')};
handles.settings.mutype = str;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
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

function edit_gaussianstd_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussianstd = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_gaussianx_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussianx = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_gaussiany_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.gaussiany = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_mhlower_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhlower = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_mhupper_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhupper = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_mhT_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.settings.mhT = val;
[mu, x1, x2] = generatemu(handles);
handles.settings.mu = mu;
handles.settings.x1 = x1;
handles.settings.x2 = x2;
rendermu(handles);
guidata(hObject, handles);

function edit_xupper_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
if val > handles.settings.xlim(1)
    set(hObject, 'String', num2str(val));
    handles.settings.xlim(2) = val;
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
    guidata(hObject, handles);
end

function edit_xlower_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
if val < handles.settings.xlim(2)
    set(hObject, 'String', num2str(val));
    handles.settings.xlim(1) = val;
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
    guidata(hObject, handles);
end

function edit_yupper_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
if val > handles.settings.ylim(1)
    set(hObject, 'String', num2str(val));
    handles.settings.ylim(2) = val;
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
    guidata(hObject, handles);
end

function edit_ylower_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
if val < handles.settings.ylim(2)
    set(hObject, 'String', num2str(val));
    handles.settings.ylim(1) = val;
    [mu, x1, x2] = generatemu(handles);
    handles.settings.mu = mu;
    handles.settings.x1 = x1;
    handles.settings.x2 = x2;
    rendermu(handles);
    guidata(hObject, handles);
end

function cbox_spherical_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.settings.spherical = val;
guidata(hObject, handles);

function cbox_datetitle_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.settings.datetitle = val;
guidata(hObject, handles);

function edit_starttime_Callback(hObject, eventdata, handles)

val = datenum(get(hObject, 'String'));
set(hObject, 'String', datestr(val, 0));
handles.settings.starttime = val;
guidata(hObject, handles);

function cbox_overwrite_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.overwrite = val;
guidata(hObject, handles);

function cbox_main_enable_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.main.enable = val;
guidata(hObject, handles);

function edit_main_rate_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.rate = val;
guidata(hObject, handles);

function edit_main_width_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.width = val;
guidata(hObject, handles);

function edit_main_height_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.main.height = val;
guidata(hObject, handles);

function edit_main_filename_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
handles.outputsettings.main.filename = str;
guidata(hObject, handles);

function cbox_main_animation_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.main.animation = val;
guidata(hObject, handles);

function cbox_cov_enable_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.coverage.enable = val;
guidata(hObject, handles);

function edit_cov_rate_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.rate = val;
guidata(hObject, handles);

function edit_cov_width_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.width = val;
guidata(hObject, handles);

function edit_cov_height_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.coverage.height = val;
guidata(hObject, handles);

function edit_cov_filename_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
handles.outputsettings.coverage.filename = str;
guidata(hObject, handles);

function cbox_cov_animation_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.coverage.animation = val;
guidata(hObject, handles);

function cbox_mu_enable_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.mu.enable = val;
guidata(hObject, handles);

function edit_mu_rate_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.rate = val;
guidata(hObject, handles);

function edit_mu_width_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.width = val;
guidata(hObject, handles);

function edit_mu_height_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mu.height = val;
guidata(hObject, handles);

function edit_mu_filename_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
handles.outputsettings.mu.filename = str;
guidata(hObject, handles);

function cbox_mu_animation_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.mu.animation = val;
guidata(hObject, handles);

function cbox_con_enable_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.convergence.enable = val;
guidata(hObject, handles);

function edit_con_rate_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.rate = val;
guidata(hObject, handles);

function edit_con_width_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.width = val;
guidata(hObject, handles);

function edit_con_height_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.convergence.height = val;
guidata(hObject, handles);

function edit_con_filename_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
handles.outputsettings.convergence.filename = str;
guidata(hObject, handles);

function cbox_con_animation_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.convergence.animation = val;
guidata(hObject, handles);

function cbox_mh_enable_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.mesohyperbolicity.enable = val;
guidata(hObject, handles);

function edit_mh_rate_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.rate = val;
guidata(hObject, handles);

function edit_mh_width_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.width = val;
guidata(hObject, handles);

function edit_mh_height_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.height = val;
guidata(hObject, handles);

function edit_mh_filename_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
handles.outputsettings.mesohyperbolicity.filename = str;
guidata(hObject, handles);

function cbox_mh_animation_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value');
handles.outputsettings.mesohyperbolicity.animation = val;
guidata(hObject, handles);

function edit_mh_period_Callback(hObject, eventdata, handles)

str = get(hObject,'String');
val = eval(str);
set(hObject, 'String', num2str(val));
handles.outputsettings.mesohyperbolicity.T = val;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UI element creation functions, all boilerplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_N_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_umax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_tstop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ntargets_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_K_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_findradius_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_findtimeconst_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ngrid_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_h_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_nsamplepts_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_lambda_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popup_algorithm_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_v1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_v2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_agentuncertainty_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_targetuncertainty_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popup_mutype_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussianstd_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussianx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gaussiany_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhlower_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhupper_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mhT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_xupper_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_xlower_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_yupper_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ylower_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_starttime_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_width_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_height_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_main_filename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_width_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_height_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_filename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_width_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_height_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_con_filename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_width_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_height_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mu_filename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_filename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_height_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_width_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cov_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mh_period_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
