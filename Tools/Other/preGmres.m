function [x,flag,relres,iter,resvec] = preGmres(fmm,b,restart,tol,maxit,x,varargin)
%GMRES   Generalized Minimum Residual Method.


if (nargin < 2)
    error(message('MATLAB:gmres:NumInputs'));
end


m = size(b,1);
n = m;
if ~iscolumn(b)
    error(message('MATLAB:gmres:Vector'));
end


% Assign default values to unspecified parameters
if (nargin < 3) || isempty(restart) || (restart == n)
    restarted = false;
else
    restarted = true;
    restart = max(restart, 0);
end
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    if restarted
        maxit = min(ceil(n/restart),10);
    else
        maxit = min(n,10);
    end
end
maxit = max(maxit, 0);

if restarted
    outer = maxit;
    if restart > n
        warning(message('MATLAB:gmres:tooManyInnerItsRestart',restart, n));
        restart = n;
    end
    inner = restart;
else
    outer = 1;
    if maxit > n
        warning(message('MATLAB:gmres:tooManyInnerItsMaxit',maxit, n));
        maxit = n;
    end
    inner = maxit;
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                   % Norm of rhs vector, b
if (n2b == 0)                    % if    rhs vector is all zeros
    x = zeros(n,1);              % then  solution is all zeros
    flag = 0;                    % a valid solution has been obtained
    relres = 0;                  % the relative residual is actually 0/0
    iter = [0 0];                % no iterations need be performed
    resvec = 0;                  % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('gmres',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 6) && ~isempty(x))
    if ~isequal(size(x),[n,1])
        error(message('MATLAB:gmres:XoSize', n));
    end
else
    x = zeros(n,1);
end

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;
minupdated = 0;

x0iszero = (norm(x) == 0);
r = b - fmm.computeVel(x);
normr = norm(r);                 % Norm of initial residual
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = [0 0];
    resvec = normr;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end
minv_b = b;

normr = norm(r);                 % norm of the preconditioned residual
n2minv_b = norm(minv_b);         % norm of the preconditioned rhs
clear minv_b;
tolb = tol * n2minv_b;
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2minv_b;
    iter = [0 0];
    resvec = n2minv_b;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

U = zeros(n,inner);
R = zeros(inner,inner);
w = zeros(inner+1,1);

for outiter = 1 : outer
    %  Construct u for Householder reflector.
    %  u = r + sign(r(1))*||r||*e1
    u = r;
    normr = norm(r);
    beta = scalarsign(r(1))*normr;
    u(1) = u(1) + beta;
    u = u / norm(u);
    
    U(:,1) = u;
    
    %  Apply Householder projection to r.
    %  w = r - 2*u*u'*r;
    w(1) = -beta;
    
    for initer = 1 : inner
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*(u(initer)')*u;
        v(initer) = v(initer) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        for k = (initer-1):-1:1
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        v = v/norm(v);
        
        %  Apply A to v.
        if initer == 1
            v = fmm.computeVel(v);
        else 
            v = fmm.reducedComputeVel(v);
        end

        %  Form Pj*Pj-1*...P1*Av.
        for k = 1:initer
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        
        %  Determine Pj+1.
        if (initer ~= length(v))
            %  Construct u for Householder reflector Pj+1.
            u = v;
            u(1:initer) = 0;
            alpha = norm(u);
            if (alpha ~= 0)
                alpha = scalarsign(v(initer+1))*alpha;
                %  u = v(initer+1:end) +
                %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
                u(initer+1) = u(initer+1) + alpha;
                u = u / norm(u);
                U(:,initer+1) = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                v(initer+2:end) = 0;
                v(initer+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:initer-1
            tmpv = v(colJ);
            v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
            v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
        end
        
        %  Compute Given's rotation Jm.
        if ~(initer==length(v))
            rho = norm(v(initer:initer+1));
            J(:,initer) = v(initer:initer+1)./rho;
            w(initer+1) = -J(2,initer).*w(initer);
            w(initer) = conj(J(1,initer)).*w(initer);
            v(initer) = rho;
            v(initer+1) = 0;
        end
        
        R(:,initer) = v(1:inner);
        
        normr = abs(w(initer+1));
        resvec((outiter-1)*inner+initer+1) = normr;
        normr_act = normr;
        
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
            if evalxm == 0
                ytmp = R(1:initer,1:initer) \ w(1:initer);
                additive = U(:,initer)*(-2*ytmp(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + ytmp(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + ytmp(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                if norm(additive) < eps*norm(x)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                xm = x + additive;
                evalxm = 1;
            elseif evalxm == 1
                addvc = [-(R(1:initer-1,1:initer-1)\R(1:initer-1,initer))*...
                    (w(initer)/R(initer,initer)); w(initer)/R(initer,initer)];
                if norm(addvc) < eps*norm(xm)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                additive = U(:,initer)*(-2*addvc(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + addvc(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + addvc(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                xm = xm + additive;
            end
            r = b - fmm.computeVel(xm);
            if norm(r) <= tol*n2b
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            end
            minv_r = r;
            
            normr_act = norm(minv_r);
            resvec((outiter-1)*inner+initer+1) = normr_act;
            
            if normr_act <= normrmin
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                xmin = xm;
                minupdated = 1;
            end
            
            if normr_act <= tolb
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:gmres:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = [outiter, initer];
                    break;
                end
            end
        end
        
        if normr_act <= normrmin
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end         % ends inner loop
    
    if isempty(initer)
        initer = 0;
    end
    
    evalxm = 0;
    
    if flag ~= 0
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        if idx > 0 % Allow case inner==0 to flow through
            y = R(1:idx,1:idx) \ w(1:idx);
            additive = U(:,idx)*(-2*y(idx)*conj(U(idx,idx)));
            additive(idx) = additive(idx) + y(idx);
            for k = idx-1 : -1 : 1
                additive(k) = additive(k) + y(k);
                additive = additive - U(:,k)*(2*(U(:,k)'*additive));
            end
            x = x + additive;
        end
        xmin = x;
        r = b - fmm.computeVel(x);
        minv_r = r;

        normr_act = norm(minv_r);
        r = minv_r;
    end
    
    if normr_act <= normrmin
        xmin = x;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if flag == 3
        break;
    end
    if normr_act <= tolb
        flag = 0;
        iter = [outiter, initer];
        break;
    end
    minupdated = 0;
end         % ends outer loop

if isempty(outiter)
    outiter = 0;
    initer = 0;
    normr_act = normrmin;
end

% returned solution is that with minimum residual
if flag == 0
    relres = normr_act / n2minv_b;
else
    x = xmin;
    iter = [imin jmin];
    relres = normr_act / n2minv_b;
end

resvec = resvec(1:max(outiter-1,0)*inner+initer+1);
if flag == 2 && initer ~= 0
    resvec(end) = [];
end

% only display a message if the output flag is not used
if nargout < 2
    if restarted
        itermsg(sprintf('gmres(%d)',restart),tol,maxit,[outiter initer],flag,iter,relres);
    else
        itermsg(sprintf('gmres'),tol,maxit,initer,flag,iter(2),relres);
    end
end

function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end
