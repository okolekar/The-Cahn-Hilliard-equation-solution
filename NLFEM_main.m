function NLFEA_main()
% p = Gaussian points
% GW = Gaussian Weights
% e = number of elements
% cm1g = global concentration at current time step
% cmg is the global cm vector initially cm = cm1g
% dt is the Time increment
% x1,x2 = elemental global co-ordinates used for Jacobi
% The first for loop commented as ....(1) fills the cm1g as cu or cl
% depending on whether the node lies below or above L/2 length only for c
% points and not for c' points.
    p = [-0.3399810435848456 0.3399810435848456 -0.861136311594053 0.861136311594053];
    GW = [0.652145154862546 0.347854845137454]; 
    e = input("Number of elements ");
    cm1g = zeros(e*2+2,1); 
    dt = input("Enter the time step ");
    x1 = 0; 
    x2 = (10/e); 
    for temp = 1: 2*e+2                             %.............(1)
        if (rem(temp,2) ~= 0)
            if (temp < (e+1))
                cm1g(temp) = 0.2;
            else
                cm1g(temp) = 0.8;
            end
        end
    end
    cmg = cm1g; %cm1g = global concentration field at current time
    %% Main Program
    % nrc = Newton Raphson Scheme itiration tracker
    % g = gauss quadrature loop itirator
    % W = Gaussian weight
    % zeta = Gaussian point
    % The second for loop itirates for time increment
    % The While loop ....(3) runs the Newton Raphson Scheme and the condition is
    % below with an if statement that breaks the while loop if the solution
    % converges.
    for t = 0:dt:86500                                 %..............(2)
        nrc = 1;
    while (1)                                     %...............(3)
        Res = zeros(e*2+2,1);
        K_t = zeros(e*2+2);
        for g = 1:4 
            if (g<2)
                W = GW(1);
            else
                W = GW(2);
            end
            zeta = p(g); 
            for j =1:e %elemental routine
                    A = Assembly(j,e);
                    cmp1e = A*cm1g;
                    cm = A*cmg;
                    [Re,Ktele] = msubroutine(zeta,x1,x2,cmp1e,cm,dt); %Fint elemental calculation Material routine
                    Re = W*Re;
                    Res = Res + A'*Re;
                    if nrc == 1
                        Fint = Res;
                    end
                    K_t = K_t + (A')*Ktele*A;
            end
        end
        delc = linsolve(K_t,Res);
        cm1g = delc + cm1g;
        nrc = nrc+1;
        if (norm(Res)< 0.005*(max(max(abs(Fint),10^-8))) && max(abs(delc)) < max(abs(cmg)) || nrc >10)
            plot(cmp1e)
            break
        end
    end
    cmg=cm1g;
    end
end

function [D2energy, D3energy,Mobility,Dmob,Ddmob] = internal(c)
    Mo = 0.3;
    psio = 0.04;
    X = 1.3;
    Mobility = 64 * Mo * c^4*(1-c)^4;
    Dmob = 256*Mo*(c^3 *(1-c)^4-c^4*(1-c)^3);
    alpha = (3*c^2*(1-c)^4-4*c^3*(1-c)^3-4*c^3*(1-c)^3+3*c^4*(1-c)^2);
    Ddmob = 265*Mo*alpha;
    D2energy = psio*X*(2-12*c+12*(c^2));
    D3energy = psio*X*(-12+24*c);
end

function As = Assembly(p,n)% n = total number of elements and p = current element
    As = zeros(4,n*2+2);
    c = 2*p-1;
    As(1,c) = 1;
    As(2,c+1) = 1;
    As(3,c+2) = 1;
    As(4,c+3) = 1;
end

function [Res,Ktelem] = msubroutine(zeta,x1,x2,cmp1,cm,delt)
    Jacobi = ((x2-x1)/2);%Jaicobi.
    he = 0.01;
    Lbd = 0.01; %Lambda
    N = [0.25*(2+zeta)*(zeta-1)^2, 0.125*he*(zeta+1)*(zeta-1)^2, 0.25*(2-zeta)*(zeta+1)^2,0.125*he*(zeta-1)*(zeta+1)^2];
    B = [0.75*(zeta^2 - 1),0.125*he*(3*zeta^2-2*zeta-1),0.75*(1-zeta^2),0.125*he*(3*zeta^2+2*zeta-1)]*(1/Jacobi);
    G = [(3/2*zeta),0.125*he*(6*zeta-2),-3/2*zeta,0.125*he*(6*zeta+2)]*(1/Jacobi)*(1/Jacobi);
    [ddpsi,d3psi,M,dM,ddM] = internal(N*cmp1);
    Q = (M*ddpsi*(B')*B*cmp1*Jacobi + dM*Lbd*G*cmp1*(B')*B*cmp1*Jacobi + M*Lbd*(G')*G*cmp1*Jacobi);
    Res = ((N')*N*(cmp1-cm)*Jacobi + delt*Q);
    D1 = (N')*N*Jacobi; 
    D2 = (M*ddpsi*(B')*B*Jacobi + M*d3psi*(B')*B*cmp1*Jacobi*N + dM*ddpsi*(B')*B*cmp1*Jacobi*N);  
    D3 = (dM*Lbd*G*cmp1*(B')*B*Jacobi + dM*Lbd*G*(B')*B*cmp1 + ddM*Lbd*G*cmp1*(B')*B*cmp1*N*Jacobi);
    D4 = (M*Lbd*(G')*G*Jacobi + dM*Lbd*(G')*G*cmp1*Jacobi*N);
    Ktelem = (D1 + D2 + D3 +D4);
end
