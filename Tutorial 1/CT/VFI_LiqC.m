%Optimized for speed by SeHyoun Ahn

clear all; clc;

tic;

s = 2; %CRRA utility with parameter s
r = 0.03; %interest rate
% s = 0.5;
% r = 0.045;
rho = 0.05; %discount rate
z1 = .1;
z2 = .2;
z = [z1,z2];
la1 = 0.02;
la2 = 0.03;
la = [la1,la2];


I=500;
amin = -0.02; %borrowing constraint
amax = 2;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
zz = ones(I,1)*z;

%b(a) and omega(a) functions
%it doesn't really matter what these are except for which a's is b(a)=0
omega = ones(I,2);
i_crit = 100;
omega(1:i_crit,:) = 0;
%bb = omega.*aa;
%bmin = 0;

fraction = 1;
liqinc = zz + fraction*r.*aa;

maxit= 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

%INITIAL GUESS
v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;

% z_ave = la2/(la1+la2)*z(1) + la1/(la1+la2)*z(2);
% v0(:,1) = (z_ave + r.*a).^(1-s)/(1-s)/rho;
% v0(:,2) = (z_ave + r.*a).^(1-s)/(1-s)/rho;

v = v0;

for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(i_crit+1:I,:) = (V(i_crit+1:I,:)-V(i_crit:I-1,:))/da;
    dVb(1:i_crit,:) = (z + r.*a(1:i_crit)).^(-s); %state constraint boundary condition
        
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
%     cf(1:i_crit,:) = min(cf(1:i_crit,:),liqinc(1:i_crit,:));
    ssf = zz + r.*aa - cf;
    Hf = cf.^(1-s)/(1-s) + dVf.*ssf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
%     cb(1:i_crit,:) = min(cb(1:i_crit,:),liqinc(1:i_crit,:));
    ssb = zz + r.*aa - cb;
    Hb = cb.^(1-s)/(1-s) + dVb.*ssb;
    %consumption and derivative of value function at steady state
    c0 = zz + r.*aa;
    dV0 = c0.^(-s);
    
    % makes a choice of forward or backward differences based on the sign of the drift    
    Ineither = (1-(ssf>0)) .* (1-(ssb<0));
    Iunique = (ssb<0).*(1-(ssf>0)) + (1-(ssb<0)).*(ssf>0);
    Iboth = (ssb<0).*(ssf>0);
    Ib = Iunique.*(ssb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(ssf>0) + Iboth.*(Hf>=Hb);
    I0 = Ineither;
    
    c = cf.*If + cb.*Ib + c0.*I0; %important to include third term
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = - Ib.*ssb/da;
    Y = - If.*ssf/da + Ib.*ssb/da;
    Z = If.*ssf/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B = (rho + 1/Delta)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

dV = dVf.*If + dVb.*Ib + dV0.*I0;

% Graphs
set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

% Verr = c.^(1-s)/(1-s) + dV_Upwind.*(zz + r.*aa - c) + ones(I,1)*la.*(V_switch - V) - rho.*V;
% 
% set(gca,'FontSize',14)
% plot(a,Verr,'LineWidth',2)
% grid
% xlabel('k')
% ylabel('Error in HJB Equation')
% xlim([amin amax])

adot = zz + r.*aa - c;

figure(1)
subplot(1,2,1)
plot(a,V,'LineWidth',2)
set(gca,'FontSize',14)
grid
xlabel('a')
ylabel('V_i(a)')
xlim([amin amax])

subplot(1,2,2)
plot(a,dV,'LineWidth',2)
set(gca,'FontSize',14)
grid
xlabel('a')
ylabel('V_i\prime(a)')
xlim([amin amax])
print -depsc value_function.eps

figure(2)
subplot(1,2,1)
plot(a,c,'LineWidth',2)
set(gca,'FontSize',14)
grid
xlabel('a')
ylabel('c_i(a)')
xlim([amin amax])

subplot(1,2,2)
plot(a,adot,a,zeros(1,I),'k--','LineWidth',2)
set(gca,'FontSize',14)
grid
xlabel('a')
ylabel('s_i(a)')
xlim([amin amax])
print -depsc policy_functions.eps