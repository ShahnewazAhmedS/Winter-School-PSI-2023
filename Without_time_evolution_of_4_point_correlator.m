clc;
clear;
close all;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PARAMETER DEFINTION FOR CALCULATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = sqrt(-1);
J = 1;
h = 1;
chi = 1;
phi= 0;

%system_sizes = 20:2:60; 
%%% uncomment the following line, if you want to computer g(U,V)
%%% in local computer please use system_size = 200 for better result

system_sizes = 100; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% INLINE FUNCTION DEFINTION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ek=@(k) 2*J*sqrt((cos(k)-h/J).^2+chi.^2 .*sin(k).*sin(k));
Zk=@(k) 2*(h-J*cos(k));
yk=@(k) 2*chi*J*sin(k);
uk=@(k) (k==0) * (1)/sqrt(2)+ (k~=0).*...
    ((Ek(k)+Zk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps);
vk=@(k) (k==0)*(-1)/sqrt(2)+(k~=0).*...
    ((i*yk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps);


gZ=@(z,zb) (sqrt(1-sqrt(1-z)).*sqrt(1-sqrt(1-zb))+...
    sqrt(1+sqrt(1-z)).*sqrt(1+sqrt(1-zb)))./(2*(1-z).^(1/8).*(1-zb).^(1/8));

gpol=@(r,theta) gZ(r*exp(i*theta),r*exp(-i*theta));


num_of_systems = length(system_sizes);
zero_sigma_sigma = zeros(1,num_of_systems);
sigma_sigma_epsilon = zeros(1,num_of_systems);
OPE_sigma_sigma_epsilon = zeros(1,num_of_systems);

% γ_k = (V_{kj} c†_j + U_{kj} c_j)
% <vac|γ_p γ_q|vac> =  Σ_{m=1)^L U_{pm} V_{qm} exp(-i*(p+q)*m) # c vacuum
% <vac|γ_p γ_q γ_r γ_s|vac> = Pfaffian((U*V')-(U*V')')
for jj = 1:num_of_systems
    
    L = system_sizes(jj);
    
    % <0|0> 
    kABC= pi/L*(kron(2*(1:(L/2))-1,[-1,1])') ;                                   % k values of Aperiodic Boundary Condition
    ABC = exp(-i*phi)/sqrt(L)*exp(-i*kABC*(1:L));                                % Fourier transformation
    uABC = conj(diag(uk(kABC)))*ABC;                                             % U matrix for all γ values, U* matrix would be V matrix for γ†
    vABC = conj(diag(vk(kABC)))*ABC;                                             % V matrix for all γ values, V* matrix would be U matrix for γ†
    mABC = triu([conj(flipud(vABC));uABC]*transpose([conj(flipud(uABC));vABC])); % [U(Π γ†) U(Π γ)]*transpose([V(Π γ†) V(Π γ)])
    mABC = mABC-transpose(mABC);                                                 % (U*V') - (U*V')'
    nABC = sqrt(real(pfaffian(mABC)));                                           % Norm of ABC vacuum
    
    % <φ_σ|φ_σ>
    kPBC = pi/L*(2*([0;kron(1:(L/2-1),[-1,1])']));                               % k values of Periodic Boundary Condition
    PBC = exp(-i*phi)/sqrt(L)*exp(-i*kPBC*(1:L));                                % Fourier transformation
    uPBC = conj(diag(uk(kPBC)))*PBC;                                             % U matrix for γ, U* matrix would be V matrix for γ†
    vPBC = conj(diag(vk(kPBC)))*PBC;                                             % V matrix for γ, V* matrix would be U matrix for γ†
    mPBC = triu([conj(flipud(vPBC));uPBC]*transpose([conj(flipud(uPBC));vPBC])); % [U(Π γ†) U(Π γ)]*transpose([V(Π γ†) V(Π γ)])
    mPBC = mPBC-transpose(mPBC);                                                 % (U*V') - (U*V')'
    nPBC = sqrt(real(pfaffian(mPBC)));                                           % Norm of PBC vacuum
    
    % <0|σ_1^x|φ_σ>
    c_ins1 = [1 zeros(1,L-1)];                                                   % σ_1^x = (V_{kj} c†_j + U_{kj} c_j) where U_{kj} = V_{kj} = 𝛿(k=j=1)
    opc1 = triu([conj(flipud(vABC));c_ins1;uPBC]*...                             % this means σ_1^x can be represented as linear combination of c†_j, c_j
        transpose([conj(flipud(uABC));c_ins1;vPBC]));                            % [U(Π γ†) U(σ_1^x) U(Π γ)]*transpose([V(Π γ†) V(σ_1^x) V(Π γ)])
    opc1 = opc1 - transpose(opc1);                                               % (U*V') - (U*V')'
    zero_sigma_sigma(jj) = abs(pfaffian(opc1)/nABC/nPBC);                        % Normalization of |φ_σ> and |0> 
    
    % <φ_ε|φ_ε>
    mEPS = triu([conj(flipud(vABC));uABC([1,2],:);conj(vABC([2,1],:));uABC]*...  % |φ_ε> = γ†_{k= -π/L} γ†_{k=π/L} |0>
        transpose([conj(flipud(uABC));vABC([1,2],:);conj(uABC([2,1],:));vABC])); % [U(Π γ†) U(γ γ) U(γ† γ†) U(Π γ)]*transpose([V(Π γ†) V(γ γ) V(γ† γ†) V(Π γ)])
    mEPS = mEPS - transpose(mEPS);
    nEPS = sqrt(real(pfaffian(mEPS)));

    % <φ_σ|σ_1^x|φ_ε>
    opc2 = triu([conj(flipud(vPBC));c_ins1;conj(vABC([2,1],:));uABC]*...
        transpose([conj(flipud(uPBC));c_ins1;conj(uABC([2,1],:));vABC]));
    opc2 = opc2 - transpose(opc2);
    sigma_sigma_epsilon(jj) = abs(pfaffian(opc2)/nEPS/nPBC);
    OPE_sigma_sigma_epsilon(jj)=sigma_sigma_epsilon(jj)/zero_sigma_sigma(jj);

    if num_of_systems == 1
        disp(strcat("Value of OPE coeficient of Primary : ",num2str(OPE_sigma_sigma_epsilon(jj),'%.8f')))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% COMPUTATION OF Δ %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recovering conformal dimension = delta
% using Least Square Approximation
% Ax = b => x = (A'A)\(A'b)

if num_of_systems > 1 
    A = [log(system_sizes)',ones(num_of_systems,1)];
    x = (A'*A)\(A'*(log(abs(zero_sigma_sigma)))');
    delta = x(1);
    figure;hold on;grid on;
    title('System size vs correlator value')
    xlabel(' L');ylabel('\langle\sigma\sigma\rangle')
    loglog(system_sizes,abs(exp(A*x)),'blue','LineWidth',3)
    scatter(system_sizes,abs(zero_sigma_sigma),'red','o','filled')    
    disp(strcat("Value of conformal dimension : ",num2str(-delta(1),'%.8f')))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% <φ_σ|σ_1^x σ_j^x|φ_σ> %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% σ^x_j = (c†_j + c_j) Π_{m=1}^{j-1} (c†_m + c_m)(c†_m - c_m)
% σ^x_j could be written as product of linear combination of c†_j, c_j
if num_of_systems == 1
    spin_op_locations = ceil(system_sizes/50):4:(ceil(49*system_sizes/50));
    num_of_op = length(spin_op_locations);
    g_eq_time = zeros(1,num_of_op);
    
    for j = 1:num_of_op
        
        jth = spin_op_locations(j);
        c_ins = kron([eye(jth-1),zeros(jth-1,L-jth+1)],[1;-1]);
        c_dag_ins = kron([eye(jth-1),zeros(jth-1,L-jth+1)],[1;1]);

        c_ins_jth = [zeros(1,jth-1), 1, zeros(1,L-jth)];
        c_dag_ins_jth = c_ins_jth;

        U = [conj(flipud(vPBC));c_ins1;c_ins_jth;c_ins;uPBC];
        V = [conj(flipud(uPBC));c_ins1;c_dag_ins_jth;c_dag_ins;vPBC];

        opc = triu(U*transpose(V));
        opc = (opc - transpose(opc));
        
        % <σσσσ> / <σσ>^2
        g_eq_time(j) = (pfaffian(opc))/(nPBC^2)/(zero_sigma_sigma(end)^2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% PLOT OF EXACT FORMULA AND ISING MODEL %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;hold on;grid on;
    title('g(u,v) at |z|=1')
    xlabel('spin operator location in angle (degree)');
    ylabel('\langle\sigma\sigma\sigma\sigma\rangle/\langle\sigma\sigma\rangle^2')
    angles = spin_op_locations/system_sizes*2*pi-2*pi/system_sizes;
    plot(linspace(min(angles),max(angles))*360/2/pi,gpol(1,linspace(min(angles),max(angles))),'blue','LineWidth',2)
    scatter(angles*360/2/pi,abs(g_eq_time),'red','o','filled')    
    axis([0,360,min(abs(g_eq_time)), max(abs(g_eq_time))])
    legend('Exact','Ising model')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% PLOT OF ERROR BETWEEN EXACT AND ISING %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;grid on;
    g_exact=gpol(1,angles);
    semilogy(angles*360/2/pi,(abs(g_exact-abs(g_eq_time))./g_exact)*100,'blue','LineWidth',2)
    xlabel('spin operator location in angle (degree)');
    ylabel('Relative Error (%) ')
    title("Error in computing g(u,v)")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PFFAFIAN FUNCTION DEFINTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pfaffian_val]= pfaffian(A)
    [n, m] = size(A);
    if n~=m
        disp("Not square matrix");
        pfaffian_val = NaN;
        return
    end
    if mod(n, 2) == 1
        pfaffian_val = 0;
        return 
    end
    pfaffian_val = 1.0;
    for k  =1:2:n
        [kmax,kind] = max(abs(A((k+1):end, k)));
        kp = k + kind;
        if kp ~= k + 1
            temp = A(k + 1, k:end);
            A(k + 1, k:end) = A(kp, k:end);
            A(kp, k:end) = temp;
            temp = A(k:end, k + 1);
            A(k:end, k + 1) = A(k:end, kp);
            A(k:end, kp) = temp;
            pfaffian_val = pfaffian_val*(-1);
        end
        if A(k + 1, k) ~= 0.0
            tau = A(k, k + 2 :end);
            tau = tau / A(k, k + 1);
            pfaffian_val = pfaffian_val* A(k, k + 1);
            if (k + 2 <= n)
                A(k + 2 :end, k + 2 :end) = A(k + 2 :end, k + 2 :end) +   transpose(tau)*transpose( A(k + 2 :end, k +1)) ; 
                A(k + 2 :end, k + 2 :end) = A(k + 2 :end, k + 2 :end) - A(k + 2 :end, k +1)*tau;
            end
        else
            pfaffian_val = 0;
        end
     end
 end
 













