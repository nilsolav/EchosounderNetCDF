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



%% Test file
% 
% ncfile='D:\repos\Github\EchosounderNetCDF\test\demo_mask.nc';
% h5disp(ncfile)
% dat0 = h5info(ncfile,'/Interpretation')
% dat0 = h5info(ncfile)
% dat = h5read(ncfile,'categories');
% 
% dat = h5read(ncfile);
% 

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
%     % write out EVL file to be used as the bottom line
%     evlFile = fullfile(evDir, [name(1:end-5) '.evl']);
%     fid = fopen(evlFile, 'w');
%     fprintf(fid, 'EVBD 3 3.00.41\r\n%d\r\n', length(bttm));
%     for j = 1:length(bttm)
%         pingDate = datestr(timestamps(j), 'yyyymmdd');
%         pingTime = datestr(timestamps(j), 'HHMMSSFFF');
%         fprintf(fid, '%s %s0 %.3f 3\r\n', pingDate, pingTime, bttm(j));
%     end
%     fclose(fid);

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
%     netcdf mask {
% 	:date_created = "20190819T134900Z";
% 	:mask_convention_version = "0.1";
% 	:mask_convention_name = "SONAR-netCDF4";
% 	:mask_convention_authority = "ICES, IMR";
% 	:rights = "Unrestricted rights";
% 	:license = "None";
% 	:Conventions = "CF-1.7, ACDD-1.3, SONAR-netCDF4-2.0";
% 	:keywords = "scrutinisation mask, echosounder";
% 	:summary = "Contains definitions of echogram scrutiny masks";
%     :title = "Echogram scrutiny masks";
	
        % write out EVR file
        fid = fopen(NCfile, 'w');
        fprintf(fid,'netcdf mask {\r\n');
        fprintf(fid,['\t\t:date_created = "',dat.date_created,'";\n']);
        fprintf(fid,['\t\t:mask_convention_version = "',dat.mask_convention_version,'";\n']);
        fprintf(fid,['\t\t:mask_convention_name = "',dat.mask_convention_name,'";\n']);
        fprintf(fid,['\t\t:mask_convention_authority = "',dat.mask_convention_authority,'";\n']);
        fprintf(fid,['\t\t:rights = "',dat.rights,'";\n']);
        fprintf(fid,['\t\t:license = "',dat.license,'";\n']);
        fprintf(fid,['\t\t:Conventions = "',dat.Conventions,'";\n']);
        fprintf(fid,['\t\t:keywords = "',dat.keywords,'";\n']);
        fprintf(fid,['\t\t:summary = "',dat.summary,'";\n']);
        fprintf(fid,['\ttitle = "',dat.title,'";\n']);
        
        % Fifth parameter is region type:
        % 0 = bad (no data), 1 = analysis, 2 = marker, 3 = fishtracks, 4 = bad (empty water)
        
        if p.Results.exportLayers
            regionId = writeRegions(fid, layer, timestamps, regionId, 1, frequency);
        end
        if p.Results.exportSchools
            regionId = writeRegions(fid, school, timestamps, regionId,  1, frequency);
        end
        if p.Results.exportErased
            regionId = writeRegions(fid, e, timestamps, regionId, 4, frequency);
        end
        if p.Results.exportExcluded
            writeRegions(fid, x, timestamps, regionId, 0, frequency);
        end
        
        fclose(fid);
    end

end


function regionId = writeRegions(fid, region, timestamps, regionId, regionType, frequency)

    regionStructureVersion = 13; % from Echoview manual
    selected = 0; % from Echoview manual
    regionCreationType = -1; % == No type
    numLinesNotes = 1;
    notes = 'Converted from LSSS scrutiny';
    numLinesDetectionSettings = 1;
    detectionSettings = 'Detected by LSSS'; % Can this be left blank?
    regionName = 'default region name';
    
    for j = 1:length(region)
        
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
        channel = NaN;
        if isfield(region(j), 'channel') && ~isempty(region(j).channel)
            if isfield(region(j).channel(1), 'frequency')
                for k = 1:length(region(j).channel)
                    if strcmp(region(j).channel(k).frequency, frequency)
                        channel = k;
                        break;
                    end
                end
            end
        end
        
        % assumes that there is only 1 species classification (if not, it
        % uses the first one).
        
        if isfield(region(j), 'channel') && ~isnan(channel)
            % if the current region has species allocated to it, use that,
            % otherwise use a 'null' species.
            if isfield(region(j).channel(channel), 'species')
                regionClassification = sprintf('{species: %s, fraction: %s, freq:%s}', ...
                    region(j).channel(channel).species(1).speciesID, ...
                    region(j).channel(channel).species(1).fraction, ...
                    region(j).channel(channel).frequency);
            else
                regionClassification = '{}';
            end
        else
            regionClassification = '{}';
        end
        
        fprintf(fid, '\r\n');
        fprintf(fid, '%d %d %d %d %d -1 0 0 0 0 0 0 0\r\n', ...
            regionStructureVersion, length(t), regionId, selected, ...
            regionCreationType);
        fprintf(fid, '%d\r\n', numLinesNotes); % allow for multiple lines
        fprintf(fid, '%s\r\n', notes);
        fprintf(fid, '%d\r\n', numLinesDetectionSettings); % allow for multiple lines
        fprintf(fid, '%s\r\n', detectionSettings);
        fprintf(fid, '%s\r\n', regionClassification);
        

        for k = 1:length(t)
            fprintf(fid, '%s %d', pointToEchoview(t(k), d(k))); % vectorise this?
        end
        fprintf(fid, '%d\r\n', regionType);
        fprintf(fid, '%s\r\n', regionName);
        
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



% 
% 
% nccreate(file,'Var1')
% 
% 
% 
% 
% keyboard
% 
% netcdf mask {
% 	:date_created = "20190819T134900Z";
% 	:mask_convention_version = "0.1";
% 	:mask_convention_name = "SONAR-netCDF4";
% 	:mask_convention_authority = "ICES, IMR";
% 	:rights = "Unrestricted rights";
% 	:license = "None";
% 	:Conventions = "CF-1.7, ACDD-1.3, SONAR-netCDF4-2.0";
% 	:keywords = "scrutinisation mask, echosounder";
% 	:summary = "Contains definitions of echogram scrutiny masks";
%     :title = "Echogram scrutiny masks";
% 	
% group: Interpretation {
% 	// SUGGESTIONS OF THINGS TO ADD:
% 	// - consider a separate implementation of layers, as per LSSS
% 	// - link to categorisation database and database version
% 	// - name of echosounder files that the data came from??
% 	
% 	types:
% 		// Note: empty_water == LSSS erased; no_data == LSSS excluded
% 		byte enum region_t {empty_water = 0, no_data = 1, analysis = 2, track = 3, marker = 4};
% 		// Storing 3D regions is not yet done, but we include the region dimension here anyway
% 		byte enum region_dim_t {twoD = 0, threeD = 1};
% 		float(*) mask_depth_t;
% 		uint64(*) mask_time_t;
% 	dimensions:
% 		regions = 3; // varies to suit data. Could also be unlimited
% 		channels = 3; // varies to suit data
% 	variables:
% 		float sound_speed;
% 			sound_speed:long_name = "Sound speed used to convert echo time into range";
% 			sound_speed:standard_name = "speed_of_sound_in_sea_water";
% 			sound_speed:units = "m/s";
% 			sound_speed:valid_min = 0.0;
% 
% 		// The bounding box of each region
% 		float min_depth(regions);
% 			min_depth:long_name = "Minimum depth for each region";
% 			min_depth:units = "m";
% 			min_depth:valid_min = 0.0;
% 		float max_depth(regions);
% 			max_depth:long_name = "Maximum depth for each regions";
% 			max_depth:units = "m";
% 			max_depth:valid_min = 0.0;
% 		uint64 start_time(regions);
% 			start_time:long_name = "Timestamp of the earliest data point in each region";
% 			start_time:units = "milliseconds since 1601-01-01 00:00:00Z";
% 			start_time:axis = "T";
% 			start_time:calendar = "gregorian";
% 			start_time:standard_name = "time";
% 		uint64 end_time(regions);
% 			end_time:long_name = "Timestamp of the latest data point in each region";
% 			end_time:units = "milliseconds since 1601-01-01 00:00:00Z";
% 			end_time:axis = "T";
% 			end_time:calendar = "gregorian";
% 			end_time:standard_name = "time";
% 			
% 		region_dim_t region_dimension; 
% 			region_dimension:long_name = "Region dimension";
% 
% 		int id(regions);
% 			id:long_name = "Dataset-unique identification number for each region";
% 		string name(regions);
% 			name:long_name = "Name of each region";
% 		string provenance(regions);
% 			provenance:long_name = "Provenance of each region"; 
% 		string comment(regions);
% 			comment:long_name = "Comment for each region";
% 		string category(regions);
% 			category:long_name = "Categorisation for each region";
% 		region_t region_type(regions);
% 			region_type:long_name = "Region type";
% 		
% 		string channel_names(channels);
% 			channel_names:long_name = "Echosounder channel names";
% 		uint region_channels(regions);
% 			region_channels:long_name = "Echosounder channels that this region applies to";
% 			region_channels:description = "Bit mask derived from channel_names (index 1 of channel_names = bit 1, index 2 = bit 2, etc). Set bits in excess of the number of channels are to be ignored.";
% 			region_channels:_FillValue = 4294967295; // 2^32-1
% 			
% 		mask_time_t mask_times(regions);
% 			mask_times:long_name = "Timestamp of each mask point";
% 			mask_times:units = "milliseconds since 1601-01-01 00:00:00Z";
% 			mask_times:axis = "T";
% 			mask_times:calendar = "gregorian";
% 			mask_times:standard_name = "time";
% 		mask_depth_t mask_depths(regions);
% 			mask_depths:long_name = "Depth pairs of mask";
% 			mask_depths:units = "m";
% 			mask_depths:valid_min = 0.0;
% 	
% 	data:
% 	    // simple example regions
% 		region_dimension = twoD;
% 		sound_speed = 1496;
% 		min_depth =  0.0, 20.5, 55.0;
% 		max_depth = 10.0, 42.0, 125.2;
% 		start_time = 13210680841000, 13210680843000, 13210680845000;
% 		end_time =   13210680847000, 13210680846000, 13210680850000;
% 		id = 1, 5, 234;
% 		name = "region1", "region2", "";
% 		provenance = "KORONA-2.6.0;LSSS", "Echoview - template ABC", "Manual inspection";
% 		comment = "", "", "whale!";
% 		category = "herring", "", "whale";
% 		region_type = analysis, empty_water, analysis;
% 		channel_names = "18kHz WBT ABC", "38kHz WBT ZYX", "120kHz GPT 123";
% 		region_channels = 5, 7, 7;
% 		mask_times = {13210680841000, 13210680842000, 13210680843000, 13210680847000}, 
% 		             {13210680843000, 13210680844000, 13210680845000, 13210680846000}, 
% 					 {13210680845000, 13210680846000, 13210680846000, 13210680847000, 13210680848000, 13210680849000, 13210680850000};
% 		mask_depths = {0.0, 5.0, 0.0, 8.0, 0.0, 10.0, 0.0, 10.0}, {20.5, 25.0, 30.5, 35.0, 35.5, 40.0, 40.0, 42.0}, {55.0, 105.0, 60.0, 80.2, 100.6, 115.0, 55.0, 107.0, 55.0, 110.0, 55.0, 115.6, 55.0, 125.2};	
% 	}
% }