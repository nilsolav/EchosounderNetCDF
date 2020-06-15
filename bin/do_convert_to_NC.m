%% Initialize

dat.date_created = '20190819T134900Z';
dat.mask_convention_version = '0.1';
dat.mask_convention_name = 'SONAR-netCDF4';
dat.mask_convention_authority = 'ICES, IMR';
dat.rights = 'Unrestricted rights';
dat.license = 'None';
dat.Conventions = 'CF-1.7, ACDD-1.3, SONAR-netCDF4-2.0';
dat.keywords = 'scrutinisation mask, echosounder';
dat.summary = 'Contains definitions of echogram scrutiny masks';
dat.title = 'Echogram scrutiny masks';


% Path to the example data
whr = which('LSSSreader_readsnapfiles');
[dr,~,~] = fileparts(whr);
dr = dr(1:end-3);

% Recursively list relevant data in the example directory. This can also be
% used to search files in any folder structure
files.snap = rdir(fullfile(dr,'exampledata','**','*.snap'));
files.work = rdir(fullfile(dr,'exampledata','**','*.work'));
files.raw = rdir(fullfile(dr,'exampledata','**','*.raw'));

% Match the corresponding snap, work and raw files (by file name)
files=LSSSreader_pairfiles(files);

pl = true; % Set to false for plotting the masks only (without background echograms)
%pl = false;

for file=1%:size(files.F,1)
    snap = files.F{file,1};
    work = files.F{file,2};
    raw  = files.F{file,3};
    % If no snap files are present, use the work files instead
    if isempty(snap)
        snap=work;
    end
    
    NCfile = [snap(1:end-4),'.h5'];
    
    % Store in working dir for now
    [~,f1,f2]=fileparts(NCfile);
    NCfile = [f1,f2];
    
    % Run conversion
    numRegions = convertWorkToNC(snap, raw, NCfile, 38000,dat);
end
