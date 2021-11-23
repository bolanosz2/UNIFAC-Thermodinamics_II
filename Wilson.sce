function W = WILSON(Mat,T,c) 
    a = [1,2,3;1,2,3;1,2,3] //Constantes de Wilson
    V = [1,2,3] // Volúmenes molares a T
    // Crear matriz para Aij
    A = zeros(c,c)
    for i=1:c
        for j=1:c
            if i=!j
                A[i][j]=V[i]/V[j]*(exp(-a[i][j]/(R*T)))
            else
                A[i][j]=1
    end
    suma1=0
    suma2=0
    suma3=0
    for i =1:c
        for k=1:c
            for j=1:c
                suma1 =+ (Mat[j][5]*A[k][j])
            end
            suma2 =+ (Mat[k][5]*A[k][i]/suma1)
        end
        for j=1:c
            suma3 =+ (Mat[j][5]*A[i][j])
        end
        Mat[i][3]=exp(1-suma3-suma2)
    end
    W = Mat
endfunction
