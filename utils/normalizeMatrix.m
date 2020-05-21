function [normalizedMatrix] = normalizeMatrix(inputMatrix)

sumOfAllElements = sum(sum(inputMatrix));
normalizedMatrix = inputMatrix / sumOfAllElements;

end