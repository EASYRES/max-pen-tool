% Number of PV units
num_pv = 15;

% Number of network buses
num_bus = 78;

% Create 50 random allocation patterns of PV units along the network
num_allct = 50;

alloc_pv = zeros(num_allct,num_pv);

for i = 1 : num_allct
    for j = 1 : num_pv
        while 1
            temp = randi(num_bus);
            if isempty(find(alloc_pv(i,1:j) == temp))
                alloc_pv(i,j) = temp + 1;
                break;
            end
        end
    end
end

save('Allocation pattern','alloc_pv','num_pv','num_allct')