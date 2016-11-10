system()

% p = 1 a single node



% p = 3
%fprintf('p = 3, b = %s\n', b)
phi = 0;
for i=1:s
    phi = phi + b(i)*c(i)*c(i);
end
simplify(phi)

phi = 0;
for i=1:s
    for j=1:s
        phi = phi + b(i)*a(i, j)*c(j);
    end
end
simplify(phi)




%fprintf('p = 3, bt = %s\n', bt)
phi = 0;
for i=1:s
    phi = phi + bt(i)*c(i)*c(i);
end
fprintf('computing for bt')
simplify(phi)

phi = 0;
for i=1:s
    for j=1:s
        phi = phi + bt(i)*a(j, i)*c(j);
    end
end

simplify(phi)

