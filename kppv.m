%--------------------------------------------------------------------------
% ENSEEIHT - 1SN - Analyse de donnees
% TP4 - Reconnaissance de chiffres manuscrits par k plus proches voisins
% fonction kppv.m
%--------------------------------------------------------------------------
function [Partition, Confusion] = kppv(DataA,DataT,labelA,labelT,K,ListeClass)

[Na,~] = size(DataA);
[Nt,~] = size(DataT);

Nt_test = Nt/10; % A changer, pouvant aller jusqu'� Nt

% Initialisation du vecteur d'etiquetage des images tests
Partition = zeros(Nt_test,1);
Confusion = zeros(length(ListeClass));

disp(['Classification des images test dans ' num2str(length(ListeClass)) ' classes'])
disp(['par la methode des ' num2str(K) ' plus proches voisins:'])

% Boucle sur les vecteurs test de l'ensemble de l'evaluation
for i = 1:Nt_test
    
    disp(['image test n' num2str(i)])

    % Calcul des distances entre les vecteurs de test 
    % et les vecteurs d'apprentissage (voisins)
    d_xi = sum((DataA - DataT(i,:)).^2,2);
    
    % On ne garde que les indices des K + proches voisins
    [~, indices] = sort(d_xi);
    indices = indices(1:K);
    
    % Comptage du nombre de voisins appartenant a chaque classe
    classes = hist(labelA(indices),ListeClass);
    
    % Recherche de la classe contenant le maximum de voisins
    
    
    % Si l'image test a le plus grand nombre de voisins dans plusieurs  
    % classes differentes, alors on lui assigne celle du voisin le + proche,
    % sinon on lui assigne l'unique classe contenant le plus de voisins 
    
    
    % Assignation de l'etiquette correspondant � la classe trouvee au point 
    % correspondant a la i-eme image test dans le vecteur "Partition" 
    [~, indice_max] = max(classes);
    Partition(i) = indice_max - 1;
    
    % Matrice de confusion
    Confusion(labelT(i)+1,Partition(i)+1) = Confusion(labelT(i)+1,Partition(i)+1) + 1;
    
end

