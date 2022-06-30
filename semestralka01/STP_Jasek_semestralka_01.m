close  all
clear all
clc

%% Příklad č. 1
%matice prechodu
U = [0,1,0,0,0,0;
     0,0,1,0,0,0;
     0,0,0,1,0,0;
     0,0,0,0,1,0;
     0,0,0,0,0,1;
     1,0,0,0,0,0]
 
U = [1/2,1/2,0,0,0,0;
     0,1/2,1/2,0,0,0;
     0,0,1/2,1/2,0,0;
     0,0,0,1/2,1/2,0;
     0,0,0,0,1/2,1/2;
     1/2,0,0,0,0,1/2]

%matice stredniho poctu kroku
syms m11 m12 m13 m14 m15 m16...
     m21 m22 m23 m24 m25 m26...
     m31 m32 m33 m34 m35 m36...
     m41 m42 m43 m44 m45 m46...
     m51 m52 m53 m54 m55 m56...
     m61 m62 m63 m64 m65 m66

M = [m11 m12 m13 m14 m15 m16;
     m21 m22 m23 m24 m25 m26;
     m31 m32 m33 m34 m35 m36;
     m41 m42 m43 m44 m45 m46;
     m51 m52 m53 m54 m55 m56;
     m61 m62 m63 m64 m65 m66]
 
M_overline = M
M_overline(1:(1+size(M,1)):end) = 0

E = ones(6)

reseni = solve(M==U*M_overline+E, [m11 m12 m13 m14 m15 m16 ...
                                   m21 m22 m23 m24 m25 m26 ...
                                   m31 m32 m33 m34 m35 m36 ...
                                   m41 m42 m43 m44 m45 m46 ...
                                   m51 m52 m53 m54 m55 m56 ...
                                   m61 m62 m63 m64 m65 m66])
                              
M = [reseni.m11 reseni.m12 reseni.m13 reseni.m14 reseni.m15 reseni.m16;
     reseni.m21 reseni.m22 reseni.m23 reseni.m24 reseni.m25 reseni.m26;
     reseni.m31 reseni.m32 reseni.m33 reseni.m34 reseni.m35 reseni.m36;
     reseni.m41 reseni.m42 reseni.m43 reseni.m44 reseni.m45 reseni.m46;
     reseni.m51 reseni.m52 reseni.m53 reseni.m54 reseni.m55 reseni.m56;
     reseni.m61 reseni.m62 reseni.m63 reseni.m64 reseni.m65 reseni.m66]
P=U;
for i=0:1:100
    P=P*U;
end
disp("Finalni matice ppsti")
P

%% Příklad č. 2
%homogenni markovsky retezec se dvema absorpcnimi stavy
U = [1,0,0,0,0,0;
     1/3,0,0,0,2/3,0;
     0,2/3,0,1/3,0,0;
     0,0,0,1,0,0;
     0,0,0,1/3,0,2/3;
     1/3,0,2/3,0,0,0]
%stredni pocet prechodu tranzientnim stavem
I = diag([1,1,1,1])
%vynechani radku a sloupcu s absorpcnimi stavy
Q = [U([2,3,5,6],[2,3,5,6])]
T = inv(I-Q)

%stredni doba stravena v tranzientnich stavech
t = T*[1,1,1,1]'

%pravdepodobnost, ze skonci v absorpcnim stavu
R=U([2,3,5,6],[1,4])%matice ppsti ze stavu mimo absorpcniho do dvou absorpcnich
d = T*R


