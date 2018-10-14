function [size_of_cell] = cell_size(input, axis)
	tmp = cell2mat(cellfun(@size,input(1:length(input)),'uni',false));
	size_of_cell = tmp(:, axis);
end
