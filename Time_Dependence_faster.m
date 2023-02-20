clc;
clear;
close all;
format long;
i = sqrt(-1);
J = 1;
h = 1;
chi = 1;
tic
phi= 0;
Ek=@(k) 2*J*sqrt((cos(k)-h/J).^2+chi.^2 .*sin(k).*sin(k));
Zk=@(k) 2*(h-J*cos(k));
yk=@(k) 2*chi*J*sin(k);
uk=@(k) (k==0) * (1)/sqrt(2)+ (k~=0).*...
    ((Ek(k)+Zk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps);
vk=@(k) (k==0)*(i)/sqrt(2)+(k~=0).*...
    ((i*yk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps);

uk2=@(k) (k==0)* 0+ (k~=0 & k~=pi).*...
    ((Ek(k)+Zk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps)+ (k==pi).* 1;
vk2=@(k) (k==0)* 1 +(k~=0 & k~=pi).*...
    ((i*yk(k)))./sqrt(2*(Ek(k)+Zk(k)).*Ek(k)+eps) + (k==pi).* 0;

gZ=@(z,zb) (sqrt(1-sqrt(1-z)).*sqrt(1-sqrt(1-zb))+...
    sqrt(1+sqrt(1-z)).*sqrt(1+sqrt(1-zb)))./(2*(1-z).^(1/8).*(1-zb).^(1/8));

gpol=@(r,theta) gZ(r*exp(i*theta),r*exp(-i*theta));

system_sizes = 160;
num_of_systems = length(system_sizes);
zero_sigma_sigma = zeros(1,num_of_systems);
sigma_sigma_epsilon = zeros(1,num_of_systems);
OPE_sigma_sigma_epsilon = zeros(1,num_of_systems);


L = system_sizes;
jj=1;
% ABC matrix <0|0>
kABC= pi/L*(kron(2*(1:(L/2))-1,[-1,1])') ;
ABC = exp(-i*phi)/sqrt(L)*exp(-i*kABC*(1:L));
uABC = conj(diag(uk(kABC)))*ABC;
vABC = conj(diag(vk(kABC)))*ABC;
mABC = triu([conj(flipud(vABC));uABC]*transpose([conj(flipud(uABC));vABC]));
mABC = mABC-transpose(mABC);
nABC = sqrt(real(pfaffian(mABC)));
% PBC matrix <φ_σ|φ_σ>
kPBC = pi/L*(2*([0;kron(1:(L/2-1),[-1,1])']));
PBC = exp(-i*phi)/sqrt(L)*exp(-i*kPBC*(1:L));
uPBC = conj(diag(uk(kPBC)))*PBC;
vPBC = conj(diag(vk(kPBC)))*PBC;
mPBC = triu([conj(flipud(vPBC));uPBC]*transpose([conj(flipud(uPBC));vPBC]));
mPBC = mPBC-transpose(mPBC);
nPBC = sqrt(real(pfaffian(mPBC)));
% <0|σ_1^x|φ_σ>
c_ins1 = [1 zeros(1,L-1)];
opc1 = triu([conj(flipud(vABC));c_ins1;uPBC]*...
    transpose([conj(flipud(uABC));c_ins1;vPBC]));
opc1 = opc1 - transpose(opc1);
zero_sigma_sigma(jj) = abs(pfaffian(opc1)/nABC/nPBC);
%%

% <φ_σ|σ_1^x σ_j^x|φ_σ>
% σ^x_j = (c†_j + c_j) Π_{m=1}^{j-1} (c†_m + c_m)(c†_m - c_m)

kPBC0 = pi/L*(2*(((-L/2+1):(L/2))'));
PBC0 = exp(-i*phi)/sqrt(L)*exp(-i*kPBC0*(1:L));

% next we want to find c_m, c†_m in term of γ_k γ†_k where k in PBC 
% [γ ; γ†] = [U V; V* U*][c ; c†]

uPBC0 = conj(diag(uk2(kPBC0)))*PBC0;
vPBC0 = conj(diag(vk2(kPBC0)))*PBC0;


BigPBC0 = [(uPBC0), (vPBC0); conj(vPBC0), conj(uPBC0)];
BigABC0 = [(uABC), (vABC); conj(vABC), conj(uABC)];
% [c ; c†] = [X Y; Y* X*][γ ; γ†] 
InvBigPBC0 = inv(BigPBC0);

xPBC0 = InvBigPBC0(1:L,1:L);           
yPBC0 = InvBigPBC0(1:L,(L+1:end));

U4sigmaX =  xPBC0 + conj(yPBC0);
V4sigmaX =  yPBC0 + conj(xPBC0);

U4sigmaY =  conj(yPBC0) - xPBC0 ; % not really sigmaY but i*sigmaY
V4sigmaY =  conj(xPBC0) - yPBC0;

U4sigmaZ = kron(U4sigmaX,[1,0]')+ kron(U4sigmaY,[0,1]');
V4sigmaZ = kron(V4sigmaX,[1,0]')+ kron(V4sigmaY,[0,1]');

% representation of γ_k γ†_k in ABC in term of γ_k γ†_k where k in PBC 
% [γABC ; γ†ABC] = [U V; V* U*][X Y; Y* X*][γPBC ; γ†PBC]
% [γABC ; γ†ABC] = [R S; S* R*][γPBC ; γ†PBC]

gInvBigPBC0 = BigABC0/BigPBC0;

rPBC0 = gInvBigPBC0(1:L,1:L);
sPBC0 = gInvBigPBC0(1:L,(L+1:end));

U4gp2 =  rPBC0 + conj(sPBC0);
V4gp2 =  sPBC0 + conj(rPBC0);

U4untry = kron(eye(L),[1,1]'); 


% spin operator placed over |z| ring as radial quantization
spin_op_locations = ceil(system_sizes/32);%ceil(system_sizes/50):4:(ceil(49*system_sizes/50));
num_of_op = length(spin_op_locations);
g_eq_time = zeros(1,num_of_op);

time_slices =  linspace(-1,1.5,20);%(log(linspace(0.005,4,20)));
num_of_times = length(time_slices);
g_eq_angle = zeros(1, num_of_times);

for j = 1:num_of_times
    
    tau = time_slices(j);
    jth = spin_op_locations;

    V4Luntry = kron(diag(exp(tau*Ek(kPBC0))),[1,0]') + kron(eye(L),[0,1]');
    V4Runtry = kron(diag(exp(-tau*Ek(kPBC0))),[1,0]') + kron(eye(L),[0,1]');

    U4gp1 =  rPBC0 + diag(exp(tau*Ek(kABC)))*conj(sPBC0);
    V4gp1 =  sPBC0 + diag(exp(tau*Ek(kABC)))*conj(rPBC0);

    U4gpL = kron(U4gp1,[1,0]') + kron(U4gp2,[0,1]');
    V4gpL = kron(V4gp1,[1,0]') + kron(V4gp2,[0,1]');
    
    U = [U4sigmaX(1,:); U4gpL;  U4sigmaZ(1:2*(jth-1),:); U4sigmaX(jth,:);  U4untry];
    V = [V4sigmaX(1,:); V4gpL;  V4sigmaZ(1:2*(jth-1),:); V4sigmaX(jth,:);  V4Runtry];

    opc = triu(U*transpose(V));
    opc = (opc - transpose(opc));
    
    % <σσσσ> / <σσ>^2
    g_eq_angle(j) = (pfaffian(opc));
    g_eq_angle(j) = real(g_eq_angle(j)/(zero_sigma_sigma(end)^2));
end

%%
figure;hold on;grid on;
title('g(u,v) at different |z|')
xlabel('time |z|');
ylabel('\langle\sigma\sigma\sigma\sigma\rangle/\langle\sigma\sigma\rangle^2')
plot((time_slices),real(g_eq_angle),'red')    



toc
















