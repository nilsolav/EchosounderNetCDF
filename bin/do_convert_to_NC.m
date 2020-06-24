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

%Path to ncgen
str0 = ['"C:\Program Files\netCDF 4.7.4\bin\ncgen" -b '];

% group: Interpretation {
%     group: v1 { // subsequent versions of this interpretation get put in new subgroups, using the numbering system v1, v2, etc.
%         // SUGGESTIONS OF THINGS TO ADD:
%         // - consider a separate implementation of layers, as per LSSS
%         // - link to categorisation database and database version
%         // - name of echosounder files that the data came from??
% 
%         :version = "1"; // increasing integers
%         :version_save_date = "20190903T154023Z"; // ISO8601 format
%         :version_author = "GJM";
%         :version_comment = "Initial scrutiny";
dat.group(1).version = '1';
dat.group(1).version_save_date = datestr(now,30);
dat.group(1).version_author = 'GJM';
dat.group(1).version_comment = 'Initial scrutiny';

dat.data(1).region_provenance = 'Converted from LSSS scrutiny';

dat.group(1).region_provenance = "LSSS";
dat.group(1).region_comment = "";

dat.group(1).region_type = "analysis";

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
files = LSSSreader_pairfiles(files);

pl = true; % Set to false for plotting the masks only (without background echograms)
%pl = false;

for file=[1:7 9:size(files.F,1)]%LSSSreader fails on file 8. Blame Gavin.
    snap = files.F{file,1};
    work = files.F{file,2};
    raw  = files.F{file,3};
    % If no snap files are present, use the work files instead
    if isempty(snap)
        snap=work;
    end
    
    NCfile = [snap(1:end-5),'.cdl'];
    disp(NCfile)
    % Create cdl file
    numRegions = convertWorkToNC(snap, raw, NCfile, 38000,dat);
    % Create NC file
    str = [str0,NCfile];
    system(str);
end
