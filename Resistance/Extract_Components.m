function [x1, x2, x3] = ExtractComponents(x)

N = numel(x);
x1 = x(1:3:N);
x2 = x(2:3:N);
x3 = x(3:3:N);
end

