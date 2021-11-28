function [x] = Interleave_Components(x1, x2, x3)

if isrow(x1)
   x = [x1;x2;x3];
   x = x(:)';
else
    x = [x1, x2, x3]';
    x = x(:);
end
end

