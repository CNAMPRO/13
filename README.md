# Systèmes linéaires - Factorisation LU
## Vecteur de départ 
```Matlab
  A = [2 1;1 3];
  b = [1;0];
```
On cherche à résoudre Ax = b
## Fonction principale

```Matlab
  function x=syslin(A,b)
    [L,U] = mylu(A);
    y=descente(L,b)';
    x=remontee(U,y)';
  end
```
On commence par récupérer les deux matrices triagulaire L et U 
  - système triangulaire inférieur Ly = b :  avec que des 1 sur la diagonale
  - système triangulaire supérieur Ux = y
  
Ensuite on va se servir de la fonction descente pour obtenir le vecteur y

Et enfin se servir de la fonction remontée pour obtenir notre vecteur x

## Mylu

```Matlab
  function [L,U]=mylu(A)
    [n,m]=size(A);
    L=eye(n);
    for k = 1:n+1
      for i = k+1:n
        L(i,k) = A(i,k)/A(k,k);
        A(i,:)=A(i,:)-L(i,k)*A(k,:);
      end
    end
    U= A(:,1:n);
  end
```

- ### Explication des lignes
  - On récupère la taille de la matrice 
  ```Matlab
    [n,m]=size(A);
  ```
  - On commence à construire le vecteur L, on lui affecte une matrice identité de taille n
  ```Matlab
    L=eye(n);
  ```
  - Pour chaque ligne 'k' de la matrice et pour chaque ligne 'i' de la matrice en dessous de k
  ```Matlab
    for k = 1:n+1
      for i = k+1:n
  ```
  - On sauvegarde la coefficient du produit(qui permet d'avoir un zero)
  ```Matlab
    L(i,k) = A(i,k)/A(k,k);
  ```
  - Pour chaque valeur de la ligne i, la valeur va être affecté par lui même moins le produit du coefficient et la valeur de la ligne k
  ```Matlab
    A(i,:)=A(i,:)-L(i,k)*A(k,:);
  ```
  - Et enfin on affecte à U la matrice A
  ```Matlab
    U= A(:,1:n);
  ```
- ### Résultat 
  - Matrice triangulaire inférieur L
  ```bash
  L =

   1.00000   0.00000
   0.50000   1.00000
  ```
  - Matrice triangulaire supérieur U
  ```bash
  U =

   2.00000   1.00000
   0.00000   2.50000
  ```
  
**Nous voilà maintenant avec nos deux matrices L et U qui vont nous permettre de résoudre ce système**

## Descente 
  ```Matlab
    function y=descente(L,b)
      y(1)=b(1)/L(1,1);
      for i=2:length(b)
       y(i)=(b(i)-dot(L(i,1:i-1),y(1:i-1)))/L(i,i);
      end
    end
  ```
- ### Explication des lignes
  - On résout y(1)
  ```Matlab
    y(1)=b(1)/L(1,1);
  ```
  - Pour toutes les lignes de la matrice sauf la première ( que l'on a dejà résolu)
  ```Matlab
    for i=2:length(b)
  ```
  - On résout chaque equation : y prend la valeur de b moins le produit scalaire du vecteur L(i,1:i-1) et du vecteur y(1:i-1) le tout sur la valeur de L à la position i,i
  ```Matlab
    y(i)=(b(i)-dot(L(i,1:i-1),y(1:i-1)))/L(i,i);
  ```
- ### Résultat 
  Matrice y 
  ```bash
  y =

   1.00000  -0.50000
  ```
  Ne pas oublier de l'inverser pour obtenir une matrice colonne  !
  ```bash
  y = y'
  y =

   1.00000
  -0.50000
  ```
  
## Remontée 

```Matlab
  function x=remontee(U,y)
  n=length(y);
  x(n)=y(n)/U(n,n);
    for i=n-1:-1:1
      x(i)=(y(i)-sum(U(i,i+1:n).*x(i+1:n)))/U(i,i);
    end
  end
```
- ### Explication des lignes 
  - on sauvegarde dans n la taille du vecteur y 
  ```Matlab
    n=length(y);
  ```
  - on résout la première équation 
  ```Matlab
    x(n)=y(n)/U(n,n);
  ```
  - Pour toutes les lignes de la matrice sauf la dernière ( que l'on a dejà résolu) 
  ```Matlab
    for i=n-1:-1:1
  ```
  - On résout chaque equation : x prend la valeur de y moins le produit scalaire du vecteur U(i,i+1:n) et du vecteur x(i+1:n) le tout sur la valeur de U à la position i,i
  ```Matlab
    x(i)=(y(i)-sum(U(i,i+1:n).*x(i+1:n)))/U(i,i);
  ```
  
  Nous voilà donc maintenant avec notre vecteur x qui est la solution du système Ax = b
- ### Résultat 
  Matrice x 
  ```bash
  x =

   0.60000  -0.20000
  ```
  Ne pas oublier de l'inverser pour obtenir une matrice colonne  !
  ```bash
  x = x'
  x =

   0.60000
  -0.20000  
  ```
  
  ## Vérification
  
  Vérifier votre résultat sur matlab 
  ```bash
    A\b
  ```
  Resultat 
  ```bash
   x = A\b
   x =

    0.60000
   -0.20000
  ```
  
# Obtenir la matrice inverse avec la factorisation LU
Le but est d'obtenir la matrice inverse de A
On rappelle A :
```matlab
  A = [2 1;1 3];
``` 
On va donc utiliser cette fonction
```matlab
  function invA=myinv(A)
    [n,m] = size(A);
    [L,U] = mylu(A);
    for j=1:n
      b=zeros(n);
      b(j)=1;
      y=descente(L,b);
      invA(:,j)=remontee(U,y);
    end
  end
```
- ## Explication des lignes
  - n et m prennent la taille de la matrice A
  ```matlab
    [n,m] = size(A);
  ```
  - on récupère les matrices L et U par la fonction mylu vu au dessus
  ```matlab
    [L,U] = mylu(A);
  ```
  - on va résoudre autant de système qu'il y a de lignes 
  ```matlab
    for j=1:n
  ```
  - Pour chaque système, on va créer une matrice b de taille n avec que des 0 et mettre un 1 à la ligne j de la matrice
  ```matlab
    b=zeros(n);
    b(j)=1;
  ```
  - on résout le système avec nos fonctions et on stocke les résultats dans la colonne j d'une matrice
  ```matlab
    y=descente(L,b);
    invA(:,j)=remontee(U,y);
  ```
- ## Résultat
  Matrice inverse de A 
  ```bash
    invA =

       0.60000  -0.20000
      -0.20000   0.40000
  ```
- ## Vérification
  Utiliser la fonction de matlab pour vérifier le résultat
  ```bash
     invA = inv(A)
     invA =

       0.60000  -0.20000
      -0.20000   0.40000
  
  ```
  
  
  
# Solution Cramer

```Matlab
function x = cramer(A,b)
%CRAMER	Solve linear system by Cramer's Rule.
%	x = CRAMER(A,b) solves the square system A*x = b.

if det(A) == 0
   error('Matrix is singular')
end
[n,n] = size(A);
for j = 1:n
   B = A;
   B(:,j) = b;
   x(j) = det(B)/det(A);
end
x = x';
```
  
  
# Bon courage ;) 
### Réalisé par RAIMON Dylan
 
  
  
 
  


