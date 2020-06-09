function dgTime = NTTime2Mlab(ntDatetime)
    %
    
    % Convert NT time (number of 100 nanoseconds since January 1 1601 into
    % Matlab time. Loses some precision!!!!
    ntSecs = double(ntDatetime) / 10000000;
    dgTime = datenum(1601, 1, 1, 0, 0, ntSecs);
end

