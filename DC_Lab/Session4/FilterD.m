function filtered_data = FilterD(Data,header)
    % a function that calcluate correlation of Data 
    % and its header
    filtered_data = upfirdn(Data,fliplr(header));
end