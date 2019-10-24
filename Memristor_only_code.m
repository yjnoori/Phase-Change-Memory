%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**Welcome to Noori's memristor array simulator**%    
%**Memristor devices arrays**%
%If you found this code useful for your work, the authors would appreciate
%your citation of the paper: Y J Noori, C H de Groot "Modelling Resistive
%and Phase Change Memory with Passive Selector Arrays- A Matlab Tool"
%arXiv:1910.05836 (2019). Doing so, will ensure that the authors continue
%on updating the code and keep it open-access.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Measurement_V] = Code_Measured_RvsI_S(I_S)
% % 
% file = matfile('Random10x10bits')
% R = file.R

R(1:100,1:100) = 0

R_h = 1000000
R_l = 10000
Read_row = 1
Read_coloumn = 1
R(Read_row, Read_coloumn) = 0
R_Sense = 1000
R_WL = 1
R_BL = 1

m = size (R,1);
n = size (R,2);

error = 2
pass = 0.0001

R_S_WL1(1:m) = 10
R_S_WL1(Read_row) = 10
R_S_BL1(1:n) = 10
R_S_BL1(Read_coloumn) = R_Sense
R_S_WL2(1:m) = 10^8
R_S_BL2(1:n) = 10^8

V_APP_WL1(1:m) = 0.5
V_APP_WL1(Read_row) = 1
V_APP_BL1(1:n) = 0.5
V_APP_BL1(Read_coloumn) = 0
V_APP_WL2(1:m) = 0
V_APP_WL2(Read_row) = 0
V_APP_BL2(1:n) = 0
V_APP_BL2(Read_coloumn) = 0

%-------------------------------------------------

parfor i = 1:m
        for j = 1:n
            V(i,j) = V_APP_WL1(i) - V_APP_BL1(j)
        end
end
ct = 0

while error > pass
    
    parfor i = 1:m
        for j = 1:n
            if R(i,j) == 0
                
                    RR(i,j) = R_h
                   
            else
                
                    RR(i,j) = R_l
                
            end            
        end
    end
    
    
    %Define A

    tv = cell(1,m);
    parfor i = 1:m
        A = sparse((zeros(n,n)))
        for j = 1:n

            if j == 1

               A(j,j) = (1/(R_S_WL1(i)))+1/(RR(i,j))+1/R_WL
               A(j,j+1) = -1/R_WL

            elseif j == n
               A(j,j) = (1/(R_S_WL2(i)))+1/(RR(i,j))+1/R_WL
               A(j,j-1) = -1/R_WL   


            else
               A(j,j-1) = -1/R_WL;
               A(j,j) = (1/RR(i,j))+(2/R_WL)
               A(j,j+1) = -1/R_WL 

            end   

        end
        tv{i} = A;
    end
    AA = blkdiag(tv{:})

    %Define B

    tv = cell(1,m);
    parfor i = 1:m
        B = sparse((zeros(n,n)))
        for j = 1:n

            B(j,j) = -1/(RR(i,j))

        end

        tv{i} = B;
    end
    BB = blkdiag(tv{:})

    %Define C

    tv = cell(1,n);
    parfor j = 1:n
        C = sparse(zeros(n,m*n))
        for i = 1:m

            C(i,n*(i-1)+j) = 1/(RR(i,j))

        end

        tv{j} = C;
    end
    CC = vertcat(tv{:})

    %Define D, might need double checking

    tv = cell(1,n); 
    parfor j = 1:n
        D = sparse((zeros(n,m*n)))
        for i = 1:m

            if i == 1

               D(i,j) = (-1/R_S_BL1(j))+(-1/R_BL)+(-1/RR(i,j))
               D(i, n*i+j) = 1/R_BL

            elseif i == m

               D(i,n*(i-2)+j) = 1/R_BL
               D(i, n*(i-1)+j) = (-1/R_S_BL2(j))+(-1/RR(i,j))+(-1/R_BL)

            else

              D(i,n*(i-2)+j) = 1/R_BL  
              D(i, n*(i-1)+j) = (-1/R_BL)+(-1/RR(i,j))+(-1/R_BL)
              D(i, n*i+j) = 1/R_BL

            end   

        end

        tv{j} = D;
    end
    DD = vertcat(tv{:})

    %Define EW

    tv = cell(1,m);
    parfor i = 1:m
        EW = sparse((zeros(n,1)))
        for j = 1:n

            if j == 1

                EW(j,1) = V_APP_WL1(i)/R_S_WL1(i)

            elseif j == n

                EW(j,1) = V_APP_WL2(i)/R_S_WL2(i)

            else
                EW(j,1) = 0
            end

        end
        tv{i} = EW;
    end
    EWEW = vertcat(tv{:})

    %Define EB

    tv = cell(1,n); 
    parfor j = 1:n
        EB = sparse(zeros(m,1))
        for i = 1:m

            if i == 1

                EB(i,1) = -V_APP_BL1(j)/R_S_BL1(j)

            elseif i == m

                EB(i,1) = -V_APP_BL2(j)/R_S_BL2(j)

            else
                EB(i,1) = 0
            end

        end
        tv{j} = EB; 
    end
    EBEB = vertcat(tv{:})


    VV = [AA   BB
         CC   DD]  \  [EWEW
                       EBEB]

    %Reshaping the voltage vector into a square matrix of bitline and word line
    %voltages

    V_WL_sq = sparse((reshape(VV(1:m*n),[m,n])).')
    V_BL_sq = sparse((reshape(VV(((m*n)+1):(2*m*n)),[m,n])).')


    %Calculate the currents throughout the array
    II = sparse(zeros(m,n))
    PP = sparse(zeros(m,n))

    parfor i = 1:m
        I = sparse((zeros(1,n)))
        P = sparse((zeros(1,n)))
        for j = 1:n



                I(j) = (V_WL_sq(i,j) - V_BL_sq(i,j))/RR(i,j)

                P(j) = I(j)*(V_WL_sq(i,j) - V_BL_sq(i,j))


        end
        II(i,:) = I
        PP(i,:) = P  
    end
    V_temp = abs((V_WL_sq - V_BL_sq)-V)
    error = max(V_temp(:))
    V = V_WL_sq - V_BL_sq
  ct = ct +1  
end

    Per_II = (II/max(II(:)))*100
    Per_PP = (PP/max(PP(:)))*100
    imagesc(PP)

    %Measured voltage at sense resistor
    Measurement_V = (sum(II(1:n, Read_coloumn))) * R_Sense
    Measurement_V = full(Measurement_V)
    
    Measured_RR = 1/(sum(II(1:n, Read_coloumn)))

