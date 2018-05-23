clear;
close all;
load donnees;
load exercice_1;
figure('Name','Image tiree aleatoirement','Position',[0.2*L,0.2*H,0.6*L,0.5*H]);

% Seuil de reconnaissance a regler convenablement
s = 10;

% Pourcentage d'information 
per = 0.95;

% Tirage aleatoire d'une image de test :
individu = randi(37);
posture = randi(6);
chemin = './Images_Projet_2018';
fichier = [chemin '/' num2str(individu+3) '-' num2str(posture) '.jpg'];
Im=importdata(fichier);
I=rgb2gray(Im);
I=im2double(I);
image_test=I(:)';
 

% Affichage de l'image de test :
colormap gray;
imagesc(I);
axis image;
axis off;

% Nombre N de composantes principales a prendre en compte 
% [dans un second temps, N peut etre calcule pour atteindre le pourcentage
% d'information avec N valeurs propres] :
N = 2;

% N premieres composantes principales des images d'apprentissage :
C = Xc*W(:,1:N);

% N premieres composantes principales de l'image de test :

Ic = image_test - individu_moyen;
Ci = Ic*W(:,1:N);

% Determination de l'image d'apprentissage la plus proche (plus proche voisin) :

label_indiv = repmat(numeros_individus,4,1);
label_indiv = label_indiv(:);
K = 5;

% Calcul des distances entre les vecteurs de test 
% et les vecteurs d'apprentissage (voisins)
d_xi = sum((C - Ci).^2,2);

% On ne garde que les indices des K + proches voisins
[~, indices] = sort(d_xi);
indices = indices(1:K);
distance_min = min(d_xi)

% Comptage du nombre de voisins appartenant a chaque classe
classes = hist(label_indiv(indices),numeros_individus);

% Recherche de la classe contenant le maximum de voisins
[~, indice_max] = max(classes);
Partition = indice_max;


% Affichage du resultat :
if distance_min < s
	individu_reconnu = numeros_individus(indice_max);
	title({['Posture numero ' num2str(posture) ' de l''individu numero ' num2str(individu)];...
		['Je reconnais l''individu numero ' num2str(individu_reconnu)]},'FontSize',20);
else
	title({['Posture numero ' num2str(posture) ' de l''individu numero ' num2str(individu)];...
		'Je ne reconnais pas cet individu !'},'FontSize',20);
end
