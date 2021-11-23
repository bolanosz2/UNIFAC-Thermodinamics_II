//(fila, columna)
t = 100 // °C
T = t + 273.15 //K
c = 3 //Numero de componentes
XY = 0.8 //Porcentaje del maximo en el componente mas volatil

//A, B, C
Antoine = [13.7819,2726.81,217.572; //benceno
13.9320,3056.96,217.625; //tolueno
15.3144,3212.43,182.739] //1-butanol

//z, Psat (t), gamma, k, x, y
values = [0.5, 0, 0, 0, 0, 0; //benceno
0.2, 0, 0, 0, 0, 0; //tolueno
0.3, 0, 0, 0, 0, 0;] //1-butanol

//Calcula las presiones de vapor
for i=1:c
    values(i,2)=exp(Antoine(i,1)-Antoine(i,2)/(t+Antoine(i,3)))
end

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

function W = WILSON(Mat,d,T,c) 
    R = 1.987 // cal/(K mol)
    a = [0,(-245.7/4.184),160.12;(316.2/4.184),0,104.68;817.67,887.80,0] //Constantes de Wilson cal/gmol
    V = [0,0,0] // Volúmenes molares a T
    
    // Interpolación de V para 1
    if T<323.15 then
        V(1)=92.263-(92.263-86.783)*((323.15-T)/(323.15-273.15))
    elseif T>=323.15 then
        V(1)=98.537-(98.537-92.263)*((373.15-T)/(373.15-323.15))
    end
    
    //Interpolación de V para 2
    if T<353.15 then
        V(2)=113.717-(113.717-107.415)*((353.15-T)/(353.25-303.15))
    elseif T>=353.15 then
        V(2)=120.879-(120.879-113.717)*((400-T)/(400-353.15))
    end
    
    //Interpolación de V para 3
    if T<343.15 then
        V(3)=97.8-(97.8-89.873)*((343.25-T)/(343.15-273.15))
    elseif T>=343.15 then
        V(3)=108.7-(108.7-97.8)*((413.15-T)/(413.15-343.15))
    end
    // Crear matriz para Aij
    A = zeros(c,c)
    for i=1:c
        for j=1:c
            if i<>j
                A(i,j)=V(i)/V(j)*(exp(-a(i,j)/(R*T)))
            else
                A(i,j)=1
            end
        end
    end
    suma1=0
    suma2=0
    suma3=0
    for i =1:c
        for k=1:c
            for j=1:c
                suma1 =+ (Mat(j,d)*A(k,j))
            end
            suma2 =+ (Mat(k,d)*A(k,i)/suma1)
        end
        for j=1:c
            suma3 =+ (Mat(j,d)*A(i,j))
        end
        Mat(i,3)=exp(1-suma3-suma2)
    end
    W = Mat
endfunction

function W = UNIFAC(Mat,colval,colres,T,c) //Calculo de los coeficientes de actividad por el metodo de UNIFAC

    //Colocar la cantidad de cada componente
    // [k, K, R, Q, vk1, vk2, vk3] 
    vk = [1, 1, 0.9011, 0.848, 0, 0, 1; //CH3
    2, 1, 0.6744, 0.540, 0, 0, 4; //CH2
    3, 1, 0.4469, 0.228, 0, 0, 0; //CH
    4, 1, 0.2195, 0.000, 0, 0, 0; //C
    10, 3, 0.5313, 0.400, 6, 5, 0; //ACH
    12, 4, 1.2663, 0.968, 0, 1, 0; //ACCH3
    13, 4, 1.0396, 0.660, 0, 0, 0; //ACCH2
    15, 5, 1.0000, 1.200, 0, 0, 1; //OH
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
    
    W = Mat
    
    //[x,r,q,Ji,Li,lngammac,lngammar,gamma]
    info = zeros(c,8)
    for i=1:c
        info(i,1)=Mat(i,colval)
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
        W(i,colres)=info(i,8) 
    end
//    disp(info, "    xi      ri      qi          Li          Ji        ln(gammai)C    ln(gammai)R    gammai")
    
endfunction

contra = values

for i=1:2
    
    
    
    if i == 1 //RAULT sin reflujo
        
        disp("RAULT")
        
        values = contra
        
        for i=1:c
            values(i,3) = 1
        end
        Pburbuja = 0
        for i=1:c //calcula la Presión burbuja
            Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
        end
        
        Procio = 0
        for i=1:c //calcula la Presión rocio
            Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
        end
        Procio = 1/Procio
        
        P = (Procio + Pburbuja)/2 //Presión de operación
        
        VectP = [0,0,0]
        
        cont = 1
        index = 0
        valuesP = values
        Vk = 0.5
        while index == 0 & cont < 100
            
            if cont > 2
                P = Pburbuja - VectP(1)*(1-XY)*(Pburbuja-Procio)/(VectP(1)-VectP(2))
            else if cont == 1
                    P = Pburbuja
                else if cont == 2
                        P = Procio
                    end
                end
            end
            
            for i=1:c //calcula K
                valuesP(i,4) = valuesP(i,3) * valuesP(i,2)/P
            end
            
            V0 = 0.5
            V = 0
            while abs(V-V0) > 0.001
                V = V0
                F = 0
                dfdv = 0
                for i=1:c
                    F = F + valuesP(i,1)*(valuesP(i,4)-1)/(1+V*(valuesP(i,4)-1))
                    dfdv = dfdv - valuesP(i,1)*((valuesP(i,4)-1)^2)/((1+V*(valuesP(i,4)-1))^2)
                end
                
                V0 = (-F/dfdv) + V
            end
            
            
            for i=1:c
                valuesP(i,5) = valuesP(i,1)/(1+V*(valuesP(i,4)-1))
                valuesP(i,6) = valuesP(i,4) * valuesP(i,5)
            end
            
            a = 0
            for i=1:c
                if abs(valuesP(i,5)-values(i,5)) < 0.0009
                    a = a + 1
                end
                if abs(valuesP(i,6)-values(i,6)) < 0.0009
                    a = a + 1
                end
            end
            if abs(Vk-V0) < 0.0009
                a = a + 1
            end
            if a == 7
                index = 1
            end
            
            //disp(cont)
            values = valuesP
             Vk = V0
             
             if cont == 1
                VectP(1) = values(1,6)
                else if cont == 2
                    VectP(2) = values(1,6)
                end
            end
             
            //disp(P,"Presion")
            cont = cont +1
        end
        //disp(Pburbuja, "Presíon burbuja")
        //disp(Procio, "Presíon rocio")
        disp(values,"   zi      Pisat       gammai    ki         xi            yi")
        disp(P, "Presión")
        disp(V, "V")
        
        
        
        
        
    else if i == 2 //Rault con 0.95 reflujo
            
            disp("RAULT con reciclo de 0.95")
            
            values = contra
            
            for i=1:c
                values(i,3) = 1
            end
            Pburbuja = 0
            for i=1:c //calcula la Presión burbuja
                Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
            end
            
            Procio = 0
            for i=1:c //calcula la Presión rocio
                Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
            end
            Procio = 1/Procio
            
            P = (Procio + Pburbuja)/2 //Presión de operación
            
            VectP = [0,0,0]
            
            Zini = values
            cont = 1
            index = 0
            valuesP = values
            Vk = 0.5
            R = 0
            Feed = 1
            Vap = 0
            
            
            while index == 0 & cont < 100
                
                Feed2 = Feed + R
                
                Pburbuja = 0
                for i=1:c //calcula la Presión burbuja
                    Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
                end
                
                Procio = 0
                for i=1:c //calcula la Presión rocio
                    Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
                end
                Procio = 1/Procio
            
                if cont > 2
                    P = Pburbuja - VectP(1)*(1-XY)*(Pburbuja-Procio)/(VectP(1)-VectP(2))
                else if cont == 1
                        P = Pburbuja
                    else if cont == 2
                            P = Procio
                        end
                    end
                end
                
                for i=1:c //calcula K y Z despues del reciclo entrante
                    valuesP(i,1) = (Zini(i,1)*Feed + R*valuesP(i,6))/Feed2
                    valuesP(i,4) = valuesP(i,3) * valuesP(i,2)/P
                end
                
                V0 = 0.5
                V = 0

                while abs(V-V0) > 0.001
                    V = V0
                    F = 0
                    dfdv = 0
                    for i=1:c
                        F = F + valuesP(i,1)*(valuesP(i,4)-1)/(1+V*(valuesP(i,4)-1))
                        dfdv = dfdv - valuesP(i,1)*((valuesP(i,4)-1)^2)/((1+V*(valuesP(i,4)-1))^2)
                    end
                    
                    V0 = (-F/dfdv) + V
                end
                
                
                for i=1:c
                    valuesP(i,5) = valuesP(i,1)/(1+V*(valuesP(i,4)-1))
                    valuesP(i,6) = valuesP(i,4) * valuesP(i,5)
                end
                
                a = 0
                for i=1:c
                    if abs(valuesP(i,5)-values(i,5)) < 0.0009
                        a = a + 1
                    end
                    if abs(valuesP(i,6)-values(i,6)) < 0.0009
                        a = a + 1
                    end
                end
                if abs(Vk-V0) < 0.0009
                    a = a + 1
                end
                if a == 7
                    index = 1
                end
                
                //disp(cont)
                values = valuesP
                 Vk = V0
                 
                 if cont == 1
                    VectP(1) = values(1,6)
                    else if cont == 2
                        VectP(2) = values(1,6)
                    end
                end
                 
                 Vap = Feed2 * V
                 R = Vap * 0.95
                 
                 
                //disp(P)
                cont = cont +1
            end
            //disp(Pburbuja, "Presíon burbuja")
            //disp(Procio, "Presíon rocio")
            disp(values,"   zi           Pisat       gammai    ki         xi            yi")
            disp(P, "Presión")
            disp(V, "V")
            
            
            
            
        else if i == 3 //Wilson con 0.5 de reflujo
                
                
                values = contra
                
                disp("WIlson con 0.5 de reflujo")
                values = WILSON(values,1,T,c)
                
                Pburbuja = 0
                for i=1:c //calcula la Presión burbuja
                    Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
                end
                
                Procio = 0
                for i=1:c //calcula la Presión rocio
                    Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
                end
                Procio = 1/Procio
                
                P = (Procio + Pburbuja)/2 //Presión de operación
                
                VectP = [0,0,0]
                
                Zini = values
                cont = 1
                index = 0
                valuesP = values
                Vk = 0.5
                R = 0
                Feed = 1
                Vap = 0
                while index == 0 & cont < 100
                    
                    Feed2 = Feed + R
                    

                    Pburbuja = 0
                    for i=1:c //calcula la Presión burbuja
                        Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
                    end
                    
                    Procio = 0
                    for i=1:c //calcula la Presión rocio
                        Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
                    end
                    Procio = 1/Procio
                
                    if cont > 2
                        P = Pburbuja - VectP(1)*(1-XY)*(Pburbuja-Procio)/(VectP(1)-VectP(2))
                    else if cont == 1
                            P = Pburbuja
                        else if cont == 2
                                P = Procio
                            end
                        end
                    end
                    
                    for i=1:c //calcula K
                        valuesP(i,1) = (Zini(i,1)*Feed + R*valuesP(i,6))/Feed2
                        valuesP(i,4) = valuesP(i,3) * valuesP(i,2)/P
                    end
                    
                    V0 = 0.5
                    V = 0
                    while abs(V-V0) > 0.001
                        V = V0
                        F = 0
                        dfdv = 0
                        for i=1:c
                            F = F + valuesP(i,1)*(valuesP(i,4)-1)/(1+V*(valuesP(i,4)-1))
                            dfdv = dfdv - valuesP(i,1)*((valuesP(i,4)-1)^2)/((1+V*(valuesP(i,4)-1))^2)
                        end
                        
                        V0 = (-F/dfdv) + V
                    end
                    
                    
                    for i=1:c
                        valuesP(i,5) = valuesP(i,1)/(1+V*(valuesP(i,4)-1))
                        valuesP(i,6) = valuesP(i,4) * valuesP(i,5)
                    end
                    
                    valuesP = WILSON(valuesP,5,T,c)
                    
                    a = 0
                    for i=1:c
                        if abs(valuesP(i,5)-values(i,5)) < 0.0009
                            a = a + 1
                        end
                        if abs(valuesP(i,6)-values(i,6)) < 0.0009
                            a = a + 1
                        end
                    end
                    if abs(Vk-V0) < 0.0009
                        a = a + 1
                    end
                    if a == 7
                        index = 1
                    end
                    
                    //disp(cont)
                    values = valuesP
                     Vk = V0
                     
                     if cont == 1
                        VectP(1) = values(1,6)
                        else if cont == 2
                            VectP(2) = values(1,6)
                        end
                    end
                     
                    Vap = Feed2 * V
                    R = Vap * 0.95
                 
                    //disp(P)
                    cont = cont +1
                end
                //disp(Pburbuja, "Presíon burbuja")
                //disp(Procio, "Presíon rocio")
                disp(values,"   zi              Pisat       gammai           ki         xi            yi")
                disp(P, "Presión")
                disp(V, "V")
                
                
                
                
            else //UNIFAC sin reflujo
                
                values = contra
                
                disp("UNIFAC sin reflujo")
                values = UNIFAC(values,1,3,T,c)
                
                Pburbuja = 0
                for i=1:c //calcula la Presión burbuja
                    Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
                end
                
                Procio = 0
                for i=1:c //calcula la Presión rocio
                    Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
                end
                Procio = 1/Procio
                
                P = 117//(Procio + Pburbuja)/2 //Presión de operación
                
                VectP = [0,0,0]
                
                cont = 1
                index = 0
                valuesP = values
                Vk = 0.5
                while index == 0 & cont < 100
                    
                    Pburbuja = 0
                    for i=1:c //calcula la Presión burbuja
                        Pburbuja = Pburbuja + values(i,1)*values(i,2)*values(i,3)
                    end
                    
                    Procio = 0
                    for i=1:c //calcula la Presión rocio
                        Procio = Procio + values(i,1)/(values(i,2)*values(i,3))
                    end
                    Procio = 1/Procio
                
                    if cont > 2
                        P = Pburbuja - VectP(1)*(1-XY)*(Pburbuja-Procio)/(VectP(1)-VectP(2))
                    else if cont == 1
                            P = Pburbuja
                        else if cont == 2
                                P = Procio
                            end
                        end
                    end
                    
                    for i=1:c //calcula K
                        valuesP(i,4) = valuesP(i,3) * valuesP(i,2)/P
                    end
                    
                    V0 = 0.5
                    V = 0
                    while abs(V-V0) > 0.001
                        V = V0
                        F = 0
                        dfdv = 0
                        for i=1:c
                            F = F + valuesP(i,1)*(valuesP(i,4)-1)/(1+V*(valuesP(i,4)-1))
                            dfdv = dfdv - valuesP(i,1)*((valuesP(i,4)-1)^2)/((1+V*(valuesP(i,4)-1))^2)
                        end
                        
                        V0 = (-F/dfdv) + V
                    end
                    
                    
                    for i=1:c
                        valuesP(i,5) = valuesP(i,1)/(1+V*(valuesP(i,4)-1))
                        valuesP(i,6) = valuesP(i,4) * valuesP(i,5)
                    end
                    
                    valuesP = UNIFAC(valuesP,5,3,T,c)
                    
                    a = 0
                    for i=1:c
                        if abs(valuesP(i,5)-values(i,5)) < 0.0009
                            a = a + 1
                        end
                        if abs(valuesP(i,6)-values(i,6)) < 0.0009
                            a = a + 1
                        end
                    end
                    if abs(Vk-V0) < 0.0009
                        a = a + 1
                    end
                    if a == 7
                        index = 1
                    end
                    
                    //disp(cont)
                    values = valuesP
                     Vk = V0
                     
                     if cont == 1
                        VectP(1) = values(1,6)
                        else if cont == 2
                            VectP(2) = values(1,6)
                        end
                    end
                     
                    //disp(P)
                    cont = cont +1
                end
                //disp(Pburbuja, "Presíon burbuja")
                //disp(Procio, "Presíon rocio")
                disp(values,"   zi      Pisat       gammai           ki         xi            yi")
                disp(P, "Presión")
                disp(V, "V")
                
            end
        end
    end
    
end
