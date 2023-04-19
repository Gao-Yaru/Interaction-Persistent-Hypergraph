function [handle] = plot_bars(intervals, dim, min_filtration_value, max_filtration_value,t)

if t == 0
    type = '';
end
if t == 1
    type = 'H hat';
end
filename = sprintf('Dim%d %s barcode', dim, type);

file_format = 'png';

line_width = 0.5;

intervals = intervals';
threshold = 1e20;
epsilon = 1e-6;

max_finite_endpoint = -threshold;
min_finite_endpoint = threshold;

left_infinite_interval_found = 0;
right_infinite_interval_found = 0;


endpoints = intervals;

num_intervals = size(endpoints, 1);


for i = 1:num_intervals
    start = endpoints(i, 1);
    finish = endpoints(i, 2);
    
    if (finish >= threshold)
        right_infinite_interval_found = 1;
    end
    
    if (start <= -threshold)
        left_infinite_interval_found = 1;
    end
    
    if (finish < threshold && finish > max_finite_endpoint)
        max_finite_endpoint = finish;
    end
    
    if (start < threshold && start > max_finite_endpoint)
        max_finite_endpoint = start;
    end
    
    if (start > -threshold && start < min_finite_endpoint)
        min_finite_endpoint = start;
    end
    
    if (finish > -threshold && finish < min_finite_endpoint)
        min_finite_endpoint = finish;
    end
end



handle = figure;
hold on;

if (exist('max_filtration_value', 'var'))
    x_max = max_filtration_value;
elseif (right_infinite_interval_found)
    x_max = max_finite_endpoint + 0.2 * (max_finite_endpoint - min_finite_endpoint);
else
    x_max = max_finite_endpoint;
end

if (exist('min_filtration_value', 'var'))
    x_min = min_filtration_value;
elseif (left_infinite_interval_found)
    x_min = min_finite_endpoint - 0.2 * (max_finite_endpoint - min_finite_endpoint);
else
    x_min = min_finite_endpoint;
end

point_width = 0.006 * (x_max - x_min);


for i = 1:num_intervals
    start = endpoints(i, 1);
    finish = endpoints(i, 2);
    y = num_intervals - i + 1;
    
    if (finish >= threshold && start <= -threshold)
        line([x_min, x_max], [y, y], 'LineWidth', line_width);
        line([x_min, x_min], [y, y], 'Marker', '<', 'LineWidth', line_width);
        line([x_max, x_max], [y, y], 'Marker', '>', 'LineWidth', line_width);
    end
    
    if (finish >= threshold && start > -threshold)
        line([start, x_max], [y, y], 'LineWidth', line_width);
        line([x_max, x_max], [y, y], 'Marker', '>', 'LineWidth', line_width);
    end
    
    if (finish < threshold && start <= -threshold)
        line([x_min, finish], [y, y], 'LineWidth', line_width);
        line([x_min, x_min], [y, y], 'Marker', '<', 'LineWidth', line_width);
    end
    
    if (finish < threshold && start > -threshold)
        if (abs(finish - start) < epsilon)
            line([start - 0.5 * point_width, finish + 0.5 * point_width], [y, y], 'LineWidth', line_width);
        else
            line([start, finish], [y, y], 'LineWidth', line_width);
        end
    end
end

axis([x_min, x_max, 0, num_intervals + 1]);


set(gca,'YTick',[]);
set(gca,'XGrid','on','YGrid','on');

% title(sprintf('Dim%d %s barcode', dim, type));

% ylabel(sprintf('Dim %d', dim));



hold off;

if (exist('filename', 'var'))
    saveas(handle, filename, file_format);
end
end
