function [view_data, label, out_label]  = gen_outliers_zoo( outlier_ratios, nview )


if 2 == exist('data/zoo.mat')
load('data/zoo.mat');
else 
file_name = 'zoo.data.txt';
format = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
delim = ', ';
num = 101;
keep_cols = 2:18;
data = load_data( file_name, format, num, delim, keep_cols );
label = data(:, end);
data = data(:, 1:end-1)';
save('data/zoo.mat', 'data', 'label');
end

[ view_data, out_label, label ] = gen_outliers_multiview( data, label, outlier_ratios, nview );

% [ view1, view2, out_label, label ] = gen_outliers( data, label, outlier_ratios, nview );
% view_data{1} = view1;
% view_data{2} = view2;

% save('Zoo_Type123_demo.mat', 'view1', 'view2', 'label', 'out_label');

end

