function result = padArrayWithZeros(inputArray, desiredLength)
    % Check the length of the input array
    currentLength = length(inputArray);
    
    % If the current length is less than the desired length
    if currentLength < desiredLength
        % Calculate the number of zeros to prepend
        numZeros = desiredLength - currentLength;
        
        % Prepend zeros to the beginning of the array
        result = [zeros(1, numZeros), inputArray];
        disp('Padding the zeros at the beginning')
    else
        % If the length is already sufficient, return the array unchanged
        result = inputArray;
        disp(['No padding, since the length of mask:', num2str(currentLength)])
    end
end
