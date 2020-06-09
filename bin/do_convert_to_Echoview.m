%%
% Convert LSSS work files into Echoview .evr and .evl files

%%%
dataDir = '../data/2017842 Saga Sea';
dataDirRaw = fullfile(dataDir, 'EK60 RAWDATA');
dataDirWork = fullfile(dataDir, 'work');
dataDirEV = fullfile(dataDir, 'echoview');
d = dir(fullfile(dataDirWork, 'L*.work'));
freqs = [1 38; 2 120];
koronaRawFiles = false;

%%%
dataDir = '../data/2018815 Juvel';
dataDirRaw = fullfile(dataDir, 'KORONA');
dataDirWork = fullfile(dataDir, 'WORK');
dataDirEV = fullfile(dataDir, 'echoview');
d = dir(fullfile(dataDirWork, 'D*.work'));
freqs = [1 38; 2 70; 3 120];
koronaRawFiles = true;

%%%
dataDir = '../data/2016001 Saga Sea';
dataDirRaw = fullfile(dataDir, 'KORONA');
dataDirWork = fullfile(dataDir, 'WORK');
dataDirEV = fullfile(dataDir, 'echoview');
d = dir(fullfile(dataDirWork, 'L*.work'));
freqs = [1 38; 2 120];
koronaRawFiles = true;

for i = 1:length(d)
    if exist(fullfile(dataDirWork, d(i).name), 'file')
        disp(['Converting ' d(i).name ' (' num2str(i) ' of ' num2str(length(d)) ')'])
        r = convertWorkToEchoview(d(i).name, dataDirWork, ...
            dataDirRaw, dataDirEV, '120', ...
            'exportSchools', true, 'exportLayers', false, ...
            'channelFrequencies', freqs, ...
            'koronaRawFiles', koronaRawFiles);
        if r == 0
            disp(' No regions written')
        end
    else
        disp(['No work file found for ' d(i).name])
    end
end

%%
% and merge these into 1 evl and 1 evr file for ease of importing into
% Echoview

% EVL files
evlFilename = fullfile(dataDirEV, 'combined.evl');
if exist(evlFilename, 'file')
    delete(evlFilename)
end

d = dir(fullfile(dataDirEV, '*.evl'));
[~, ind] = sort({d.name});
d = d(ind);

j = 1;
clear l
for i = 1:length(d)
    %disp([num2str(i) ' of ' num2str(length(d))])
    fid = fopen(fullfile(d(i).folder, d(i).name), 'r');
    fgetl(fid);
    num = fscanf(fid, '%d\r\n', 1);
    for i = 1:num
        l{j} = fgetl(fid);
        j = j + 1;
    end
    fclose(fid);
end

fid = fopen(evlFilename, 'w');
fprintf(fid, '%s\r\n', 'EVBD 3 3.00.41');
fprintf(fid, '%d\r\n', length(l));
for i = 1:length(l)
    fprintf(fid, '%s\r\n', l{i});
end
fclose(fid);
disp(['Merged ' num2str(length(d)) ' .evl files'])

% EVR files
evrFilename = fullfile(dataDirEV, 'combined.evr');
if exist(evrFilename, 'file')
    delete(evrFilename)
end

d = dir(fullfile(dataDirEV, '*.evr'));
[~, ind] = sort({d.name});
d = d(ind);

j = 1;
num = 0;
clear l
for i = 1:length(d)
    %disp([num2str(i) ' of ' num2str(length(d))])
    fid = fopen(fullfile(d(i).folder, d(i).name), 'r');
    fgetl(fid);
    n = fscanf(fid, '%d', 1);
    fgetl(fid);
    num = num + n;
    while ~feof(fid)
        l{j} = fgetl(fid);
        j = j + 1;
    end
    fclose(fid);
end

fid = fopen(evrFilename, 'w');
fprintf(fid, '%s\r\n', 'EVRG 6 3.00.41');
fprintf(fid, '%d\r\n', num);
for i = 1:length(l)
    fprintf(fid, '%s\r\n', l{i});
end
fclose(fid);
disp(['Merged ' num2str(length(d)) ' .evr files'])
