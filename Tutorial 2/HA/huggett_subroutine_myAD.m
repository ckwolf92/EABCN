function objective = huggett_subroutine_myAD(r_t, data)

global N I amax amin da s z zz aa dt rho

v_st = data{1};
gg0 = data{2};
Aswitch = data{3};

SS = r_t(1) * zeros(N,1);
dS = r_t(1) * zeros(N,1);

V = v_st;
v = r_t(1) .* ones(size(V,1),size(V,2),N);
    
    for n=N:-1:1
        dVf = r_t(1) .* V;
        dVb = r_t(1) .* V;
        
        v(:,:,n)=V;
        % forward difference
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVf(I,:) = (z + r_t(n).*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
        % backward difference
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVb(1,:) = (z + r_t(n).*amin).^(-s); %state constraint boundary condition
        
        I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
        
        %consumption and savings with forward difference
        cf = dVf.^(-1/s);
        ssf = zz + r_t(n).*aa - cf;
        %consumption and savings with backward difference
        cb = dVb.^(-1/s);
        ssb = zz + r_t(n).*aa - cb;
        %consumption and derivative of value function at steady state
        c0 = zz + r_t(n).*aa;
        dV0 = c0.^(-s);
        
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift
        If = ssf > 0; %positive drift --> forward difference
        Ib = ssb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state
        %make sure backward difference is used at amax
        %Ib(I,:) = 1; If(I,:) = 0;
        %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
        %already taken care of automatically
        
        dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
        c = dV_Upwind.^(-1/s);
        u = c.^(1-s)/(1-s);
        
        %CONSTRUCT MATRIX
        X = - min(ssb,0)/da;
        Y = - max(ssf,0)/da + min(ssb,0)/da;
        Z = max(ssf,0)/da;
        
        A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
        
        %%Note the syntax for the cell array
        A_t{n} = A;
        B = (1/dt + rho)*speye(2*I) - A;
        
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
        
        V = [V_stacked(1:I),V_stacked(I+1:2*I)];
        ss = zz + r_t(n).*aa - c;
    end
    
    %plot(a,v(:,:,1),a,v(:,:,N))
    
    gg{1}=gg0;
    for n=1:N
        AT=A_t{n}';
        %Implicit method in Updating Distribution.
        gg{n+1}= (speye(2*I) - AT*dt)\gg{n};
        %gg{n+1}=gg{n}+AT*gg{n}*dt; %This is the explicit method.
        %check(n) = gg(:,n)'*ones(2*I,1)*da;
        SS(n) = gg{n}'*aa(:)*da;
        temp = gg{n+1}';
        dS(n) = temp*aa(:)*da - gg{n}'*aa(:)*da;
    end
    
    objective = dS;