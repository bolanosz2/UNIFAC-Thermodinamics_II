//Sistema Sistema (1)Fenol (2)Benceno (3)Tolueno
//Temperatura 70C
E = 0.0001 //Tolerancia
T = 70 //Temperatura en C (Antoine)
c = 3 //Numero de componentes del sistema
z = [0.2,0.3,0.5] //[z1,z2,z3] composiciones globales
x = [0,0,0] //Composicion en fase liquida
y = [0,0,0] //Composicion en fase gaseosa
gama = [1,1,1] 

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
    
function x = getnum(M,f,c) //Para matrices con 2 titulos, como en amk
    [row,col] = size(M)
    a = 0
    b = 0
    for i=1:row
        if M(i,1) == f
            a = i
        end
    end
    for  i=1:col
        if M(1,i) == c 
            b = i
        end
    end
    x = M(a,b)    
endfunction

function x = getvalc(M,f,c) //Para matrices con titulo en las filas
    [row,col] = size(M)
    for i=1:row
        if M(i,1) == f
            x = M(i,c)
        end
    end
endfunction

function x = getvalf(M,f,c) //Para matrices con titulo en las columnas
    [row,col] = size(M)
    for i=1:col
        if M(1,i) == c
            x = M(f,i)
        end
    end
endfunction

function W = UNIFAC(vec,T,c) //Calculo de los coeficientes de actividad por el metodo de UNIFAC

    //Colocar la cantidad de cada componente
    // [k, K, R, Q, vk1, vk2, vk3] 
    vk = [1, 1, 0.9011, 0.848, 0, 0, 0; //CH3
    2, 1, 0.6744, 0.540, 0, 0, 0; //CH2
    3, 1, 0.4469, 0.228, 0, 0, 0; //CH
    4, 1, 0.2195, 0.000, 0, 0, 0; //C
    10, 3, 0.5313, 0.400, 5, 6, 5; //ACH
    12, 4, 1.2663, 0.968, 0, 0, 1; //ACCH3
    13, 4, 1.0396, 0.660, 0, 0, 0; //ACCH2
    15, 5, 1.0000, 1.200, 1, 0, 0; //OH
    17, 7, 0.9200, 1.400, 0, 0, 0; //H20
    19, 9, 1.6724, 1.488, 0, 0, 0; //CH3CO
    20, 9, 1.4457, 1.180, 0, 0, 0; //CH2CO
    25, 13, 1.1450, 1.088, 0, 0, 0; //CH3O
    26, 13, 0.9183, 0.780, 0, 0, 0; //CH2O
    27, 13, 0.6908, 0.468, 0, 0, 0; //CHO
    32, 15, 1.4337, 1.244, 0, 0, 0; //CH3NH
    33, 15, 1.2070, 0.936, 0, 0, 0; //CH2NH
    34, 15, 0.9795, 0.624, 0, 0, 0; //CHNH
    41, 19, 1.8701, 1.724, 0, 0, 0; //CH3CN
    42, 19, 1.6434, 1.416, 0, 0, 0] //CH2CN
    [rowvk,colvk]=size(vk)
    
    evl = [0, 1, 3, 4, 5, 7, 9, 13, 15, 19;
    1, 0.00, 61.13, 76.5, 968.5, 1318, 476.40, 251.5, 255.7, 597;
    3, -11.12, 0.00, 167, 636.1, 903.8, 25.77, 32.14, 122.8, 212.5;
    4, -69.7, -146.8, 0.00, 803.2, 5695, -52.1, 213.1, -49.29, 6096;
    5, 156.4, 89.6, 25.82, 0.00, 353.5, 84, 28.06, 42.7, 6.712;
    7, 300, 362.3, 377.6, -229.1, 0.00, -195.4, 540.5, 168, 112.6;
    9, 26.76, 140.1, 365.8, 164.5, 472.5, 0.00, -103.6, -174.2, 481.7;
    13, 83.36, 52.13, 65.69, 237.7, -314.7, 191.1, 0.00, 251.5, -18.51;
    15, 65.33, -22.31, 223, -150, -448.2, 394.6, -56.08, 0.00, 147.1;
    19, 24.82, -22.97, -138.4, 185.4, 242.8, -287.5, 38.81, -108.5, 0.00;]
    
    W = vec
    
    //[x,r,q,Ji,Li,lngammac,lngammar,gamma]
    info = zeros(c,8)
    for i=1:c
        info(i,1)=vec(i)
    end
    
    for i = 2:c+1 //calcula ri y qi
        r = 0
        q = 0
        for j = 1:rowvk
            if vk(j,i+3) <> 0
                r = r + (vk(j,i+3)*vk(j,3))
                q = q + (vk(j,i+3)*vk(j,4))
            end
        end
        info(i-1,2) = r
        info(i-1,3) = q
    end
    
    for i=1:c //calcula Ji y Li
        valorahi = 0
        valoralla = 0
        for j=1:c
            valorahi = valorahi + info(j,1)*info(j,2)
            valoralla = valoralla + info(j,1)*info(j,3)
        end
        info(i,4) = info(i,2)/valorahi
        info(i,5) = info(i,3)/valoralla
    end
    
    for i=1:c //calcula ln(gamma)C
        info(i,6)=1-info(i,4)+log(info(i,4))-5*info(i,3)*(1-(info(i,4)/info(i,5))+log(info(i,4)/info(i,5)))
    end
    
    eki = zeros(1,c+1)
    valorahi=0
    for i=1:rowvk //crea eki
        if vk(i,5)<>0 | vk(i,6)<>0 | vk(i,7)<>0
            valorahi=valorahi+1
            eki(valorahi,1)=vk(i,1)
        end
    end
    bik = eki' //crea bik
    [roweki,coleki]=size(eki)
    for i=2:c+1 //completa eki
        for j=1:roweki
            eki(j,i)=getvalc(vk,eki(j,1),i+3)*getvalc(vk,eki(j,1),4)/info(i-1,3)
        end
    end
    
    amkTao = zeros(1,1)
    valorahi=0
    valorahi=0
    for i=1:rowvk //crea amk y tao
        if vk(i,5)<>0 | vk(i,6)<>0 | vk(i,7)<>0
            valorahi=valorahi+1
            amkTao(1+valorahi,1)=vk(i,1)
            amkTao(1,1+valorahi)=vk(i,1)
        end
    end
    [rowamkTao,colamkTao]=size(amkTao)
    
    amk=amkTao
    for i=2:colamkTao //completa amk
        for j=2:rowamkTao
            amk(j,i)=getnum(evl,getvalc(vk,amk(j,1),2),getvalc(vk,amk(1,i),2))
        end
    end
    
    tao=amkTao
    for i=2:colamkTao //completa tao
        for j=2:rowamkTao
            tao(j,i)=exp(-amk(j,i)/T)
        end
    end
    
    for j=2:c+1 //completa bij
        for i=1:roweki 
            for k=1:roweki
                m = bik(1,k)
                bik(j,i) = bik(j,i) + getvalc(eki,m,j)*getnum(tao,m,bik(1,i))
            end
        end
    end
    
    citak = zeros(2,1)
    valorahi=0
    for i=1:rowvk //crea citak
        if vk(i,5)<>0 | vk(i,6)<>0 | vk(i,7)<>0
            valorahi=valorahi+1
            citak(1,valorahi)=vk(i,1)
        end
    end
    sk = citak //crea sk
    for g=1:roweki //completa citak
        valorahi = 0
        valoralla = 0
        for i=1:c
            valorahi = valorahi + info(i,1)*info(i,3)*getvalc(eki,citak(1,g),i+1)
            valoralla = valoralla + info(i,1)*info(i,3)
        end
        citak(2,g)=valorahi/valoralla
    end
    
    for g=1:roweki //completa sk
        valorahi = 0
        for i=1:roweki
            valorahi = valorahi + citak(2,i)*getnum(tao,citak(1,i),sk(1,g))
        end
        sk(2,g)=valorahi
    end
    
    for i=1:c //calcula ln(gamma)R
        valorahi = 0
        for k=1:roweki
            valorahi = valorahi + citak(2,k)*(bik(i+1,k)/sk(2,k))-eki(k,i+1)*log(bik(i+1,k)/sk(2,k))
        end
        info(i,7) = info(i,3)*(1-valorahi)
    end
    
    for i=1:c //calcula gamma
        info(i,8) = exp(info(i,6)+info(i,7))
    end
    
    for i=1:c
        W(i)=info(i,8) 
    end
    disp(info, "    xi      ri      qi          Li          Ji        ln(gammai)C    ln(gammai)R    gammai")
    
endfunction



antoine = Psat(T,c,CA)

j = 0
index = 0
P = 0
V=0
gama = UNIFAC(z,T,c)
while index == 0 & j<30
    V1 = 0.5
    
    if j > 1
        gama = UNIFAC(x,T,c)
    end
    
    Pburb = bublP(z,antoine,gama)
    Procio = dewP(z,antoine,gama)
    Pprom = (Pburb+Procio)/2
    if j ==0
        P = Pburb + 0.001
    else if j==1
            P = (Procio+Pburb)/2
            ymax = y(2)
         else
            yiter = y(2)
            P = Pburb-((Pburb-P)*ymax*(1-0.8)/(ymax-yiter))
        end
    end
    
    //Aqui va el calculo de gamma con UNIFAC
    k = vecK(gama,P,antoine)
    
    while abs(V1-V) > E
        V =V1
        F = 0
        dF = 0
        for i=1:c
            F = F + z(i)*(k(i)-1)/(1+V*(k(i)-1))
            dF = dF - z(i)*((k(i)-1)^2)/((1+V*(k(i)-1))^2)
        end
        V1 = (-F/dF) + V
    end
    
    for i=1:c
        x(i)=z(i)/(1+V1*(k(i)-1))
        y(i)=x(i)*k(i)
    end
    
    
    
    V = V1
    j = j + 1
end
    disp("UNIFAC sin reciclo")
    disp(P,"   Presión")
    disp(x, "   x1           x2           x3")
    disp(y, "   y1           y2           y3")
    disp(V, "   Fracción Vaporizada V")











