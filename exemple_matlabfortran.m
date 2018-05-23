% création d'une matrice symétique
A = gallery('poisson', 10);
A = full(A + A');

% appel à la routine fortran subspace_iter1
% help fortran_subspace_iter_ev pour son utilisation
[V w res it] = fortran_subspace_iter_ev(A, 50, 1, 0.2, 1e-16, 2000);

% qualité des couples propres
for i=1:size(V,2)
    norm(A*V(:,i) - V(:,i)*w(i))
end
