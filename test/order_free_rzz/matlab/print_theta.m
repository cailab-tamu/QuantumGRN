function s = print_theta(A, names)
%PRINT_THETA  Format an n x n theta matrix with gene labels.
    nl = char(10);
    s = sprintf('%8s', ''); for j = 1:numel(names), s = [s sprintf('%8s', names{j})]; end
    s = [s nl];
    for i = 1:numel(names)
        s = [s sprintf('%8s', names{i}) sprintf('%8.3f', A(i,:)) nl];
    end
end
