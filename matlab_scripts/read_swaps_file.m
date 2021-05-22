%Reads .dat files resulting from Chaste simulations.
function M = read_swaps_file(filename)
fId = fopen(filename);
line = fgetl(fId);
row_counter = 1;
while ischar(line)
   row = sscanf(line,'%f');
   M(row_counter,1:length(row)) = row';
   row_counter = row_counter+1;
   line = fgetl(fId);
end

end
