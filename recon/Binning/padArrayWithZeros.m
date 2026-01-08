function result = padArrayWithZeros(inputArray, desiredLength)
    % Check the length of the input array
    currentLength = length(inputArray);
    
    % If the current length is less than the desired length
    if currentLength < desiredLength
        % Calculate the number of zeros to prepend
        numZeros = desiredLength - currentLength;
        
        % Prepend zeros to the end of the array
        result = [inputArray, zeros(1, numZeros)];
        disp('Padding the zeros at the end')
    else
        % If the length is already sufficient, return the array unchanged
        result = inputArray;
        disp(['No padding, since the length of mask:', num2str(currentLength)])
    end
end
