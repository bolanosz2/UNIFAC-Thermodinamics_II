
E = 0.00001 //Tolerancia
T = 70 //Temperatura en C (Antoine)
c = 3 //Numero de componentes del sistema
z = [0.2,0.3,0.5] //[z1,z2,z3] composiciones globales
x = [0,0,0] //Composicion en fase liquida
y = [0,0,0] //Composicion en fase gaseosa
gama = [1,1,1] 
r=0.95 //Razon de reciclo

CA= [14.4384,3507.80,175.400; 
13.7819,2726.81,217.572;
13.9320,3056.96,217.625] // Matriz con las constantes de Antoine. Fila [A,B,C] Columna [Componente i]


//Calculo de presiones de saturación
function psat = Psat(T,c,CA)
    psat = [0,0,0];
    for i=1:c
        psat(i) = exp(CA(i,1)-(CA(i,2)/(T+CA(i,3))))
    end
endfunction

//Calculo de la presion en el puunto de burbuja donde z = x
function Pb = bublP(x,antoine,gama)
    Pb = 0
    for i=1:c
        Pb = Pb + x(i)*antoine(i)*gama(i)
    end
endfunction

//Calculo de la presion en el punto de rocio donde z = y
function Pb = dewP(y,antoine,gmma)
    Pb = 0
    for i=1:c
        Pb = Pb + y(i)/(antoine(i)*gama(i))
    end
    Pb = 1/Pb
endfunction

//Calculo de los Ki
function K = vecK(gama,P,antoine)
    i = 0
    K = [0,0,0]
    for i=1:c
        K(i)= (gama(i)*antoine(i)/P)
    end
endfunction
    

antoine = Psat(T,c,CA)

disp(antoine)
j = 0
index = 0
P = 0
V=0
in = 1
zini = z
R = 0

while index == 0 & j<30
    V1 = 0.8
    in2 = in + R
    
    Pburb = bublP(z,antoine,gama)
    Procio = dewP(z,antoine,gama)
    Pprom = (Pburb+Procio)/2
    if j ==0 //Composicion 'y' maxima del componente mas volatil, se guarda en siguiente iteracion 
        P = Pburb + 0.0001
    else if j==1 
            P = (Procio+Pburb)/2
            ymax = y(2)
         else
            yiter = y(2)
            P = Pburb-((Pburb-P)*ymax*(1-0.8)/(ymax-yiter))
        end
    end
    
    for i=1:c
        z(i) = ((zini(i)*in)+R*y(i))/in2
    end
    
    //Aqui va calculo de gamma con Wilson
    k = vecK(gama,P,antoine)
    
    while abs(V1-V) > E
        V = V1
        F = 0
        dF = 0
        for i=1:c
            F = F + z(i)*(k(i)-1)/(1+V*(k(i)-1))
            dF = dF - (z(i)*((k(i)-1)^2)/((1+V*(k(i)-1))^2))
        end
        V1 = (-F/dF) + V
    end
    
    for i=1:c
        x(i)=z(i)/(1+V1*(k(i)-1))
        y(i)=x(i)*k(i)
    end
    
    V = V1
    R = in2*V*r
    j = j + 1
end
    
    disp("Raoult con reciclo")
    disp(P,"   Presión")
    disp(x, "   x1           x2           x3")
    disp(y, "   y1           y2           y3")
    disp(V, "   Fracción Vaporizada V")










