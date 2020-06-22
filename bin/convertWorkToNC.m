function numRegions = convertWorkToNC(work, raw, NCfile, frequency, dat,varargin)
% Write the matlab representation of school boxes to the NC files

    % Converts:
    % - the lower integration line in LSSS work files into sonar-NetCDF4
    % - regions into Sonar-NetCDf4
    
    p = inputParser;
    addOptional(p, 'exportExcluded', true);
    addOptional(p, 'exportErased', true);
    addOptional(p, 'exportLayers', true);
    addOptional(p, 'exportSchools', true);
    addOptional(p, 'channelFrequencies', []); % of the form [1 18; 2 38; 3 70]. Used for the erase regions.
    parse(p, varargin{:});

    % Get scrutiny
    [school, layer, exclude, erased, info] = LSSSreader_readsnapfiles(work);


    % Get timestamps of all pings in current file
    rawName =raw;
    [rawDir,rawName,ex]=fileparts(raw);
    rawName = [rawName,ex];
    timestamps = getPingTimes(rawDir,rawName);
    
    % check that number of timestamps is consistent with number of 
    % pings used in scrutiny file.
    if length(timestamps) ~= info.numberOfPings
        warning(['Raw files has ' num2str(length(timestamps)) ...
            ', while work file has ' num2str(info.numberOfPings) ' pings.'])
        
        % try to fix it
        if info.numberOfPings-1 == length(timestamps) && length(timestamps) > 2
            timestamps(end+1) = timestamps(end) + timestamps(end) - timestamps(end-1);
        end
    end
    
    % do a sanity check on the timestamps. Sometimes the timestamp is at
    % the start of time (year 1601) - usually in the first ping in a file.
    % Fix this by interpolation.
    if timestamps(1) < datenum(1980,1,1,1,1,1)
        timestamps(1) = timestamps(2) - diff(timestamps(2:3));
    end
    
    % it looks like the lsss reader gets layer .x coordinates from 0 to
    % number of pings. Seems like it should be 1 to number of pings?
    % NOTE: looks ok for schools...

    % Find the highest ping number over all layers
    maxx = 0;
    for j = 1:length(layer)
        maxx = max([maxx max(layer(j).x)]);
    end

    % Find the deepest layer boundary for each ping to use as the bottom line
    bttm = zeros(1, maxx);
    for j = 1:length(layer)
        for k = 1:length(layer(j).x)
            bttm(layer(j).x(k)) = max([bttm(layer(j).x(k)) layer(j).y(k)]);
        end
    end

    % Find and convert the layers/schools. Erased and excluded areas also
    % get converted to Echoview regions.

    % Work files give a channel number for erased regions, not the
    % frequency. We instead require the user to tell us the mapping between
    % channel numbers and frequency. This is used to work out which erase
    % masks to use.
    channelID = [];
    if ~isempty(p.Results.channelFrequencies)
        i = find(str2double(frequency) == p.Results.channelFrequencies(:,2));
        if length(i) == 1
            channelID = p.Results.channelFrequencies(i,1);
        end
    end
    e = convertMaskToPolygon(erased, channelID);
    x = convertExcludeToPolygon(exclude, timestamps);

    numRegions = 0;
    if p.Results.exportExcluded
        numRegions = numRegions + length(x);
    end
    if p.Results.exportErased
        numRegions = numRegions + length(e);
    end
    if p.Results.exportLayers
        numRegions = numRegions + length(layer);
    end
    if p.Results.exportSchools
        numRegions = numRegions + length(school);
    end

    regionId = 0;

    if numRegions > 0
        % Create and open NC file
        fid = fopen(NCfile, 'w');
        
        % Add metadata
        
        % netcdf mask {
        %     :date_created = "20190819T134900Z";
        %     :mask_convention_version = "0.1";
        %     :mask_convention_name = "SONAR-netCDF4";
        %     :mask_convention_authority = "ICES, IMR";
        %     :rights = "Unrestricted rights";
        %     :license = "None";
        %     :Conventions = "CF-1.7, ACDD-1.3, SONAR-netCDF4-2.0";
        %     :keywords = "scrutinisation mask, echosounder";
        %     :summary = "Contains definitions of echogram scrutiny masks";
        %     :title = "Echogram scrutiny masks";
	
        fprintf(fid,'netcdf mask {\r\n');
        fprintf(fid,['\t:date_created = "',dat.date_created,'";\n']);
        fprintf(fid,['\t:mask_convention_version = "',dat.mask_convention_version,'";\n']);
        fprintf(fid,['\t:mask_convention_name = "',dat.mask_convention_name,'";\n']);
        fprintf(fid,['\t:mask_convention_authority = "',dat.mask_convention_authority,'";\n']);
        fprintf(fid,['\t:rights = "',dat.rights,'";\n']);
        fprintf(fid,['\t:license = "',dat.license,'";\n']);
        fprintf(fid,['\t:Conventions = "',dat.Conventions,'";\n']);
        fprintf(fid,['\t:keywords = "',dat.keywords,'";\n']);
        fprintf(fid,['\t:summary = "',dat.summary,'";\n']);
        fprintf(fid,['\t:title = "',dat.title,'";\n\n']);
        
        % Add Layer
        
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
        fprintf(fid,'group: Interpretation {\n');
        fprintf(fid,['\tgroup: v',dat.group(1).version,'{\n']);
        fprintf(fid,['\t\t:version = "',dat.group(1).version,'";\n']);
        fprintf(fid,['\t\t:version_save_date = "',dat.group(1).version_save_date,'";\n']);
        fprintf(fid,['\t\t:version_author = "',dat.group(1).version_author,'";\n']);
        fprintf(fid,['\t\t:version_comment = "',dat.group(1).version_comment,'";\n']);

        % types:
        %     // Note: empty_water == LSSS erased; no_data == LSSS excluded
        %     byte enum region_t {empty_water = 0, no_data = 1, analysis = 2, track = 3, marker = 4};
        %     // Storing 3D regions is not yet done, but we include the region dimension here anyway
        %     byte enum region_dim_t {twoD = 0, threeD = 1};
        %     float(*) mask_depth_t;
        %     mask_depth_t(*) mask_depths_t;
        %     uint64(*) mask_time_t; // ragged array for region ping times
        
        fprintf(fid,'\t\ttypes:\n');
        fprintf(fid,'\t\t\tbyte enum region_t {empty_water = 0, no_data = 1, analysis = 2, track = 3, marker = 4};\n')
        fprintf(fid,'\t\t\tbyte enum region_dim_t {twoD = 0, threeD = 1};\n');
        fprintf(fid,'\t\t\tfloat(*) mask_depth_t;\n');
        fprintf(fid,'\t\t\tmask_depth_t(*) mask_depths_t;\n')
        fprintf(fid,'\t\t\tuint64(*) mask_time_t;\n');

        % dimensions:
        [dimensions, categories] = getDimensions(layer, school, e ,x );
        dimensions.regions = 3; %// varies to suit data. Could also be unlimited
        dimensions.channels = 3; %// varies to suit data
        dimensions.categories = 5; %// varies to suit data.
        fprintf(fid,'\t\tdimensions:\n');
        fprintf(fid,['\t\t\tregions = ',num2str(dimensions.regions),';\n']);
        fprintf(fid,['\t\t\tchannels = ',num2str(dimensions.channels),';\n']);
        fprintf(fid,['\t\t\tcategories = ',num2str(dimensions.categories),';\n']);
        
        %variables
		fprintf(fid,'\t\tvariables:\n');
        fprintf(fid,'			float sound_speed;\n');
        fprintf(fid,'				sound_speed:long_name = "Sound speed used to convert echo time into range";\n');
        fprintf(fid,'				sound_speed:standard_name = "speed_of_sound_in_sea_water";\n');
        fprintf(fid,'				sound_speed:units = "m/s";\n');
        fprintf(fid,'	\t\t\tsound_speed:valid_min = 0.0f;\n');
        fprintf(fid,'\n');
        fprintf(fid,'\t\t\t// The bounding box of each region\n');
        fprintf(fid,'\t\t\tfloat min_depth(regions);\n');
        fprintf(fid,'\t\t\t\tmin_depth:long_name = "Minimum depth for each region";\n');
        fprintf(fid,'\t\t\t\tmin_depth:units = "m";\n');
        fprintf(fid,'\t\t\t\tmin_depth:valid_min = 0.0f;\n');
        fprintf(fid,'\t\t\tfloat max_depth(regions);\n');
        fprintf(fid,'\t\t\t\tmax_depth:long_name = "Maximum depth for each regions";\n');
        fprintf(fid,'\t\t\t\tmax_depth:units = "m";\n');
        fprintf(fid,'\t\t\t\tmax_depth:valid_min = 0.0f;\n');
        fprintf(fid,'\t\t\tuint64 start_time(regions);\n');
        fprintf(fid,'\t\t\t\tstart_time:long_name = "Timestamp of the earliest data point in each region";\n');
        fprintf(fid,'\t\t\t\tstart_time:units = "milliseconds since 1601-01-01 00:00:00Z";\n');
        fprintf(fid,'\t\t\t\tstart_time:axis = "T";\n');
        fprintf(fid,'\t\t\t\tstart_time:calendar = "gregorian";\n');
        fprintf(fid,'\t\t\t\tstart_time:standard_name = "time";\n');
        fprintf(fid,'\t\t\tuint64 end_time(regions);\n');
        fprintf(fid,'\t\t\t\tend_time:long_name = "Timestamp of the latest data point in each region";\n');
        fprintf(fid,'\t\t\t\tend_time:units = "milliseconds since 1601-01-01 00:00:00Z";\n');
        fprintf(fid,'\t\t\t\tend_time:axis = "T";\n');
        fprintf(fid,'\t\t\t\tend_time:calendar = "gregorian";\n');
        fprintf(fid,'\t\t\t\tend_time:standard_name = "time";\n');
        fprintf(fid,'\t\t\t\t\n');
        fprintf(fid,'\t\t\tregion_dim_t region_dimension; \n');
        fprintf(fid,'\t\t\t\tregion_dimension:long_name = "Region dimension";\n');
        fprintf(fid,'\n');
        fprintf(fid,'\t\t\tint region_id(regions);\n');
        fprintf(fid,'\t\t\t\tregion_id:long_name = "Dataset-unique identification number for each region";\n');
        fprintf(fid,'\t\t\tstring region_name(regions);\n');
        fprintf(fid,'\t\t\t\tregion_name:long_name = "Name of each region";\n');
        fprintf(fid,'\t\t\t\tregion_name:_Encoding = "utf-8";\n');
        fprintf(fid,'\t\t\tstring region_provenance(regions);\n');
        fprintf(fid,'\t\t\t\tregion_provenance:long_name = "Provenance of each region"; \n');
        fprintf(fid,'\t\t\t\tregion_provenance:_Encoding = "utf-8";\n');
        fprintf(fid,'\t\t\tstring region_comment(regions);\n');
        fprintf(fid,'\t\t\t\tregion_comment:long_name = "Comment for each region";\n');
        fprintf(fid,'\t\t\t\tregion_comment:_Encoding = "utf-8";\n');
        fprintf(fid,'\t\t\tint region_order(regions);\n');
        fprintf(fid,'\t\t\t\tregion_order:long_name = "The stacking order of the region";\n');
        fprintf(fid,'\t\t\t\tregion_order:comment = "Regions of the same order cannot overlap";\n');
        fprintf(fid,'\t\t\tregion_t region_type(regions);\n');
        fprintf(fid,'\t\t\t\tregion_type:long_name = "Region type";\n');
        fprintf(fid,'\t\t\t\n');
        fprintf(fid,'\t\t\t// The acosutic categories. Each layer may have several categories and proportions.\n');
        fprintf(fid,'\t\t\tstring region_category_names(categories);\n');
        fprintf(fid,'\t\t\t\tregion_category_names:long_name = "Categorisation name";\n');
        fprintf(fid,'\t\t\t\tregion_category_names:_Encoding = "utf-8";\n');
        fprintf(fid,'\t\t\tfloat region_category_proportions(categories);\n');
        fprintf(fid,'\t\t\t\tregion_category_proportions:long_name = "Proportion of backscatter for the categorisation";\n');
        fprintf(fid,'\t\t\t\tregion_category_proportions:value_range = 0.0f, 1.0f;\n');
        fprintf(fid,'\t\t\tint region_category_ids(categories);\n');
        fprintf(fid,'\t\t\t\tregion_category_ids:long_name = "region_id of this categorisation and proportion";\n');
        fprintf(fid,'\t\t\t\n');
        fprintf(fid,'\t\t\tstring channel_names(channels);\n');
        fprintf(fid,'\t\t\t\tchannel_names:long_name = "Echosounder channel names";\n');
        fprintf(fid,'\t\t\t\tchannel_names:_Encoding = "utf-8";\n');
        fprintf(fid,'\t\t\tuint region_channels(regions);\n');
        fprintf(fid,'\t\t\t\tregion_channels:long_name = "Echosounder channels that this region applies to";\n');
        fprintf(fid,'\t\t\t\tregion_channels:description = "Bit mask derived from channel_names (index 1 of channel_names = bit 1, index 2 = bit 2, etc). Set bits in excess of the number of channels are to be ignored.";\n');
        fprintf(fid,'\t\t\t\tregion_channels:_FillValue = 4294967295; // 2^32-1\n');
        fprintf(fid,'\t\t\t\t\n');
        fprintf(fid,'\t\t\tmask_time_t mask_times(regions);\n');
        fprintf(fid,'\t\t\t\tmask_times:long_name = "Timestamp of each mask point";\n');
        fprintf(fid,'\t\t\t\tmask_times:units = "milliseconds since 1601-01-01 00:00:00Z";\n');
        fprintf(fid,'\t\t\t\tmask_times:axis = "T";\n');
        fprintf(fid,'\t\t\t\tmask_times:calendar = "gregorian";\n');
        fprintf(fid,'\t\t\t\tmask_times:standard_name = "time";\n');
        fprintf(fid,'\t\t\tmask_depths_t mask_depths(regions);\n');
        fprintf(fid,'\t\t\t\tmask_depths:long_name = "Depth pairs of mask";\n');
        fprintf(fid,'\t\t\t\tmask_depths:units = "m";\n');
        fprintf(fid,'\t\t\t\tmask_depths:valid_min = 0.0f;\n');
        
        % Extract region data
        
%         if p.Results.exportLayers
%             regionId = writeRegions(fid, layer, timestamps, regionId, 1, frequency);
%         end
         if p.Results.exportSchools
             data = writeRegions(fid, school, timestamps, regionId,  1, frequency, dat);
         end
         
%         if p.Results.exportErased
%             regionId = writeRegions(fid, e, timestamps, regionId, 4, frequency);
%         end
%         if p.Results.exportExcluded
%             writeRegions(fid, x, timestamps, regionId, 0, frequency);
%         end
        
        % Concatenate region data
        
        % Write region data to cdl file
        fprintf(fid,'\n\t\tdata:\n');
        fprintf(fid,'\t\t\tregion_dimension = twoD;\n');
        fprintf(fid,'\t\t\tsound_speed = 1496;\n');
        str = ['\t\t\tmin_depth =  ',num2str(data.min_depth,'%.1f, ')];
        fprintf(fid,[str(1:end-1),';\n']);
        str = ['\t\t\tmax_depth =  ',num2str(data.max_depth,'%.1f, ')];
        fprintf(fid,[str(1:end-1),';\n']);
        str = ['\t\t\tstart_time = ',num2str(data.start_time,'%i, ')];
        fprintf(fid,[str(1:end-1),';\n']);
        str = ['\t\t\tend_time = ',num2str(data.end_time,'%i, ')];
        fprintf(fid,[str(1:end-1),';\n']);
        str = ['\t\t\tregion_id = ',num2str(data.region_id,'%i, ')];
        fprintf(fid,[str(1:end-1),';\n']);
        str0 = char(strjoin(data.region_name,'","'));
        str = ['region_name = "',str0,'"'];
        fprintf(fid,'\t\t\t%s;\n',str);
        str0 = char(strjoin(data.region_provenance,'", "'));
        str = ['region_provenance = "',str0,'"'];
        fprintf(fid,'\t\t\t%s;\n',str);
        str0 = char(strjoin(data.region_comment,'", "'));
        str = ['region_comment = "',str0,'"'];
        fprintf(fid,'\t\t\t%s;\n',str);
        str0 = char(strjoin(data.region_category_names,'", "'));
        str = ['region_category_names = "',str0,'"'];
        fprintf(fid,'\t\t\t%s;\n',str);
        str = ['region_category_proportions = ',num2str(data.region_category_proportions,'%.1f, ')];
        str(end) = ';';
        fprintf(fid,'\t\t\t%s\n',str);
        str = ['region_category_ids = ',num2str(data.region_category_ids,'%i, ')];
        str(end) = ';';
        fprintf(fid,'\t\t\t%s\n',str);
        str0 = char(strjoin(data.region_type,', '));
        str = ['region_type = ',str0];
        fprintf(fid,'\t\t\t%s;\n',str);
        str0 = char(strjoin(data.channel_names,'", "'));
        str = ['channel_names = "',str0,'"'];
        fprintf(fid,'\t\t\t%s;\n',str);
        str = ['region_channels = ',num2str(data.region_channels,'%i, ')];
        str(end) = ';';
        fprintf(fid,'\t\t\t%s\n',str);
        % Iterate mask times
        for i = 1:length(data.mask_times)
            if i==1
                fprintf(fid,'\t\t\t%s','mask_times = {');
            else
                fprintf(fid,'\t\t\t%s','             {');
            end
            for j = 1:length(data.mask_times{i})
                fprintf(fid,num2str(data.mask_times{i}{j}));
                if j<length(data.mask_times{i})
                    fprintf(fid,', ');
                end
            end
            if i<length(data.mask_times)
                fprintf(fid,'},\n');
            else
                fprintf(fid,'};\n');
            end
        end
        % Iterate mask depths
        fprintf(fid,'\t\t\t%s','mask_depths = {');
        for i = 1:length(data.mask_depths)
            for j = 1:length(data.mask_depths{i})
               fprintf(fid,'{');
                str = num2str(cell2mat(data.mask_depths{i}{j}),'%.1f, ');
                fprintf(fid,str(1:end-1));
                if j<length(data.mask_times{i})
                    fprintf(fid,'}, ');
                else
                    fprintf(fid,'}');
                end
            end
            if i<length(data.mask_depths)
                fprintf(fid,'}, {');
            else
                fprintf(fid,'};');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,'\t\t}\n');
        fprintf(fid,'\t}\n');
        fprintf(fid,'}\n');
        fclose(fid);
    end
end

%mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);


function data = writeRegions(fid, region, timestamps, regionId, regionType, frequency,dat)

% test data:
data.region_dimension = 'twoD';
data.sound_speed = 1496;
data.min_depth =  [0.0, 20.5, 55.0];
data.max_depth = [10.0, 42.0, 125.2];
data.start_time = [13189164120001, 13189164121000, 13189164124000];
data.end_time =  [13189164123004, 13189164124000, 13189164131000];
data.region_id = [1, 5, 234];
data.region_name = ["region1", "region2", ""];
data.region_provenance = ["KORONA-2.6.0;LSSS", "Echoview - template ABC", "Manual inspection"];
data.region_comment = ["", "", "whale!"];
data.region_category_names = ["herring", "krill", "seal", "lion", "platypus"];
data.region_category_proportions = [0.9, 0.1, 0.45, 0.40, 0.10];
data.region_category_ids = [1, 1, 234, 234, 234];
data.region_type =  ["analysis", "empty_water", "analysis"];
data.channel_names = ["18kHz WBT ABC","38kHz WBT ZYX", "120kHz GPT 123"];
data.region_channels = [5, 7, 7];
data.mask_times = {{13189164120001, 13189164121002, 13189164122003, 13189164123004}, ...
    {13189164121000, 13189164122000, 13189164123000, 13189164124000}, ...
    {13189164124000, 13189164125000, 13189164126000, 13189164127000, 13189164128000, 13189164129000, 13189164131000}};
data.mask_depths = {{{0.0, 15.0}, {0.0, 4.0, 5.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}}, {{20.5, 25.0}, {30.5, 35.0}, {35.5, 40.0}, {40.0, 42.0}}, {{55.0, 105.0}, {60.0, 80.2, 100.6, 115.0}, {55.0, 107.0}, {55.0, 110.0}, {55.0, 115.6}, {55.0, 125.2}, {60, 115}}};
data2 = data;
clear data
% Real data

% Per region
%     region_dimension = 'twoD';
%     sound_speed = 1496;
%     min_depth = [0.0, 20.5, 55.0];
%     max_depth = [10.0, 42.0, 125.2];
%     start_time = 13189164120001;
%     end_time =  13189164123004,
%     region_id = 1;
%     region_name = "region1";
%     region_provenance = 'KORONA-2.6.0;LSSS'
%     region_comment = "";
%     region_category_names = "herring";
%     region_category_proportions = 0.9;
%     region_category_ids = 1;
%     region_type =  'analysis';
%     channel_names = '38kHz WBT ABC';
%     region_channels = 5;
%     
%     mask_times = {13189164120001, 13189164121002, 13189164122003, 13189164123004};
%     mask_depths = {{0.0, 15.0}, {0.0, 4.0, 5.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}};

% Create id variable that increment by one at each new cat/channel
category_ids =1;

% Loop over regions
for j = 1:length(region)
    numLinesNotes = 1;
    
    % get the ping interval around the start and end of the
    % region (may need it later)
    leftX = min(region(j).x);
    if leftX == 1
        leftX = 2;
    end
    
    % sometimes things are 1 ping too long, so trim
    %region(j).x(region(j).x > length(timestamps)) = length(timestamps);
    
    rightX = max(region(j).x);
    if rightX == length(timestamps)
        rightX = rightX - 1;
    end
    
    % if the ping number is larger than the ping timestamps, trim the
    % region
    region(j).x(region(j).x > length(timestamps)) = length(timestamps);
    
    % if the file has 1 ping, just simulate things a little
    if leftX <= 1 || leftX > length(timestamps)
        startPingInt = 0.1/86400; % [days]
    else
        startPingInt = diff(timestamps(leftX-1:leftX)); % [days]
    end
    if rightX <= 1 || rightX+1 > length(timestamps)
        endPingInt = 0.1/86400; % [days]
    else
        endPingInt =  diff(timestamps(rightX:rightX+1)); % [days]
    end
    
    t = timestamps(region(j).x);
    d = region(j).y;
    regionId = regionId + 1;
    
    % Fix up some undesired region stuff
    if length(d) == 4 && length(unique(region(j).x)) <= 2 && length(unique(d)) == 2
        % For rectangular regions, make them look a little better in
        % Echoview by expanding the start and stop pings to cover all
        % of that ping (otherwise they start/stop in the middle of
        % pings.
        t = [min(t)-startPingInt/2 min(t)-startPingInt/2 max(t)+endPingInt/2 max(t)+endPingInt/2];
        d = [min(d) max(d) max(d) min(d)];
    end
    
    % Work out which channel has the requested frequency. If it is not
    % there, then skip that region in the EVR file.
    
    % Channel to freq information is only available in the layer and
    % school structures. Exclude regions apply across all frequencies.
    % Erased regions are channel specific.
    
%     channel = NaN;
%     if isfield(region(j), 'channel') && ~isempty(region(j).channel)
%         if isfield(region(j).channel(1), 'frequency')
%             for k = 1:length(region(j).channel)
%                 if strcmp(region(j).channel(k).frequency, frequency)
%                     channel = k;
%                     break;
%                 end
%             end
%         end
%     end

    % assumes that there is only 1 species classification (if not, it
    % uses the first one). NOH: I assume this is a special case for
    % ev files? We need to expand and use the category_ids, right?
    
    if isfield(region(j), 'channel')
        for channel=1:length(region(j).channel)
            % if the current region has species allocated to it, use that,
            % otherwise use a 'null' species ("").
            if isfield(region(j).channel(channel), 'species')
                for sp = 1:length(region(j).channel(channel).species)
                    data.region_category_names(category_ids) = ...
                        string(region(j).channel(channel).species(sp).speciesID);
                    data.region_category_proportions(category_ids) = ...
                        str2double(region(j).channel(channel).species(sp).fraction);
                    data.region_category_ids(category_ids) = category_ids;
                    data.channel_names(category_ids) = string(region(j).channel(channel).frequency);
                    data.region_channels(category_ids) = j;
                    category_ids = category_ids + 1;
                end
            else
                data.region_category_names(category_ids) = "0";
                data.region_category_proportions(category_ids) = 1;
                data.region_category_ids(category_ids) = category_ids;
                data.channel_names(category_ids) = string(region(j).channel(channel).frequency);
                data.region_channels(category_ids) = j;
                category_ids = category_ids + 1;
            end
        end
    end    
    % Bounding box
    data.min_depth(j) = min(region(j).y);
    data.max_depth(j) = max(region(j).y);
    data.start_time(j) = interp1(1:length(timestamps),timestamps,min(region(j).x));
    data.end_time(j) = interp1(1:length(timestamps),timestamps,max(region(j).x));
    % IDdata
    data.region_id(j) = j;
    data.region_name(j) = string(['Layer',num2str(j)]);
    data.region_provenance(j) = dat.group(1).region_provenance;
    data.region_comment(j) = dat.group(1).region_comment;
    data.region_type(j) = dat.group(1).region_type;
    % Mask times
    c = unique(region(j).x);
    data.mask_times{j} = num2cell(c);
    
    % Mask depths
    
    for n=1:length(c)
        % Sort the depths within this mask time
        ind = region(j).x == c(n);
        data.mask_depths{j}{n} = num2cell(sort(region(j).y(ind)));
    end
end
end

function p = convertExcludeToPolygon(exclude, timestamps)
    % converts LSSS's exlcude region form into polygons
    
    if isempty(exclude)
        p = struct([]);
        return
    end

    
    tol = 0.2/86400; % 0.2 sec
    
    p = [];
    for i = 1:length(exclude)
        % find the ping timestamps that corresponds to the start of the
        % exclude region
        j = find(abs(timestamps - exclude(i).startTime) < tol);
        j = j(1); % just take the first if there are multiple
        
        if isempty(j)
            error(['No ping time found that matches with start of exclude region ' num2str(i)])
        end
        p(i).x = [j j j+exclude(i).numOfPings-1 j+exclude(i).numOfPings-1];
        p(i).y = [0 9999 9999 0]; % depths in metres
    end    
end

function e = convertMaskToPolygon(erased, channelID)
    % converts regions in mask format into polygons.
    
    if isempty(erased)
        e = struct([]);
        return
    end
    
    % find the largest ping in all of the masks
    max_ping = 0;
    for ch_i = 1:length(erased.channel)
        max_ping = max([max_ping max(erased.channel(ch_i).x)]);
    end
    % make up the mask from the erased data
    resolution = 0.1; % [m]
    maxRange = 1500; % [m]
    D = false(max_ping, maxRange / resolution);

    % LSSS has separate masks for each channel. If we are given a
    % channelID, just output that, otherwise merge all of the LSSS channel
    % masks together.
    
    for ch_i = 1:length(erased.channel)
        if isempty(channelID) || erased.channel(ch_i).channelID == channelID
            x = erased.channel(ch_i).x;
            y = erased.channel(ch_i).y;
        
            for i = 1:length(x)
                ranges = cumsum(y{i});
                for j = 1:2:length(ranges)
                    start = floor(ranges(j)/resolution);
                    stop  =  ceil(ranges(j+1)/resolution);
                    if start < 1
                        start = 1;
                    end
                    if stop > maxRange / resolution
                        stop = maxRange / resolution;
                    end
                
                    D(x(i), start:stop) = true;
                end
            end
        end
    end
    
    e = mask2poly(D');
    
    % ignore holes for the moment - should really split surrounding
    % polygons into two and include the holes (which won't be holes
    % anymore).
    j = 1;
    for i = 1:length(e)
        % mask2poly produces a triangular region when a mask is one ping
        % wide, so detect and fix that.
       if e(i).Length == 3 && range(e(i).X) == 0
           e(i).X(end+1) = e(i).X(end);
           e(i).Y = [min(e(i).Y) max(e(i).Y) max(e(i).Y) min(e(i).Y)];
           e(i).Length = 4;
       end
       % and we get a region with only two points with a one-ping-wide
       % region that is at the end of the data
       if e(i).Length == 2 && range(e(i).X) == 0
          e(i).X = ones(1,4) * e(1).X(1);  
          e(i).Y = [min(e(i).Y) max(e(i).Y) max(e(i).Y) min(e(i).Y)];
          e(i).Length = 4;
       end
       
       %if ee(i).isFilled % not a hole
       %    e(j) = ee(i);
       %    j = j + 1;
       %end
    end
    
    % and make polygons in the form that the rest of the code expects
    for i = 1:length(e)
        e(i).x = e(i).X;
        e(i).y = e(i).Y * resolution;
    end
    e = rmfield(e, {'X', 'Y'});
end

function s = pointToEchoview(timestamp, depth)
    % convert timestamp and depth into the Echoview .evr file form

    d = datestr(timestamp, 'yyyymmdd');
    t = datestr(timestamp, 'HHMMSSFFF');

    s = sprintf('%s %s0 %f', d, t, depth);

end

function t = getPingTimes(rawDir, fname)

    % go through the raw file and pick out the timestamps of the RAW
    % datagrams.

    headerLength = 12;
    i = 1;

    fid = fopen(fullfile(rawDir, fname));

    while(1)
        dglength = fread(fid,1,'int32');
        if feof(fid)
            break
        end
        
        type = char(fread(fid,4,'char')');    
        lowdatetime = fread(fid,1,'uint32');
        highdatetime = fread(fid,1,'uint32');
        fseek(fid, dglength-headerLength, 0); % skip the datagram data
        fread(fid, 1, 'int32'); % the trailing datagram marker
        dt = NTTime2Mlab(highdatetime*2^32 + lowdatetime);

        if ~isempty(type) && strcmp(type(1:3), 'RAW') 
            t(i) = dt;
            i = i + 1;
        end
    end
    
    fclose(fid);
    
    % Gets us every unique ping
    t = unique(t);

end

function [dim,categories] = getDimensions(layer, school, e ,x )
%dimensions(1)

%     regions = 3; // varies to suit data. Could also be unlimited
%     channels = 3; // varies to suit data
%     categories = 5; // varies to suit data.

% Sum the number of regions
dim.regions = length(layer) + length(school) + length(e) + length(x);

% Get the number of unique categories
k=1;
type = {'layer','school','e','x'};
mch=1;
for i=1:4
    eval(['region=',type{i},';'])
    for j=1:length(region)
        if isfield(region(j), 'channel') 
            % if the current region has species allocated to it, use that,
            % otherwise use a 'null' species.

            % Find the number of channels
            mch = max(mch,length(region(j).channel));
            for ch=1:length(region(j).channel)
            if isfield(region(j).channel(ch), 'species')
                cat{k}=region(j).channel(ch).species(1).speciesID;
                k=k+1;
            end
            end
        end
    end
end

dim.channels = mch;
categories = unique(cat);
dim.categories = length(categories);

end
