%%
% Reads in and plots the region definitions in the demo mask netcdf file.

f ='demo_mask.nc';

% read in some data
t = h5read(f, '/Interpretation/mask_times');
d = h5read(f, '/Interpretation/mask_depths');
mind = h5read(f, '/Interpretation/min_depth');
maxd = h5read(f, '/Interpretation/max_depth');
st = h5read(f, '/Interpretation/start_time');
et = h5read(f, '/Interpretation/end_time');
id = h5read(f, '/Interpretation/id');
name = h5read(f, '/Interpretation/name');

c = h5read(f, '/Interpretation/sound_speed');

% and some attributes
time_units = h5readatt(f, '/Interpretation/mask_times', 'units');
depth_units = h5readatt(f, '/Interpretation/mask_depths', 'units');
c_units = h5readatt(f, '/Interpretation/sound_speed', 'units');

clf
% plot each bounding box and the ping masks inside the boxes
for i = 1:length(t)
    plot([st(i) et(i) et(i) st(i) st(i)],[mind(i) mind(i) maxd(i) maxd(i) mind(i)])
    hold on
    
    text(double(st(i)), double(mind(i)), [num2str(id(i)) ', (' name{i} ')'])

    
    dd = reshape(d{i}, 2, length(d{i})/2)';
    tt = t{i};
    
    for j = 1:length(tt)
        plot([tt(j) tt(j)], [dd(j,1) dd(j,2)], 'LineWidth', 2)
    end
end

title(['Regions (c=' num2str(c) ' [' c_units '])'])

set(gca, 'Ydir', 'reverse')
xlabel(['Time (' time_units ')'])
ylabel(['Depth (' depth_units ')'])

