function scft_AB_2D()
N   = 50; %NxN box with each cell of length and breadth dr
NP  = 40; %no. of segments in diblock chain
fA  = 0.5; %fraction A in diblock
fB  = 1.0 - fA;
NA  = NP*fA;       %number of A beads
NB  = NP*fB;       %number of B beads

beta= 0.1; 
NG  = NP*beta; %NG = no. of beads in grafted chain
alpha=0.1; %volume ration of particle to diblock chain
sigma=3.0; %no. of polymer chain grafted on each particle

phigpb = 0.15; %phigpb = volume fraction of diblock copolymer
phigb  = (sigma*beta/(alpha + sigma*beta))*phigpb;
phipb  = phigpb - phigb;
phicb  = 1.0 - phigpb; %phicb = volume fraction of polymer grafted particle

sel = 0.250;  %epsilon
ChiAB = 0.13; %chiAB = incompatibility between A and B
ChiAP = sel*ChiAB;
ChiBP = (1-sel)*ChiAB;

ds  = 0.01;     %bin of chain, diblock chain
dsg = 0.01;     %bin of chain, grafted chain
NPF = (1.0/ds) + 1; % total no. of ds points
NPG = (1.0/dsg) + 1;
NSFA= (1.0/ds)*fA; %no. of points of A segment
NSFB= (1.0/ds)*fB;
dr  = 0.1;      %bin of box, N X N area
V   = dr*dr*dr*N*N*N; % area of box

lambda = 0.1; %???????????
eps    = 0.1;%????????????

for j=1:N
    for i=1:N
        for k=1:N
        w(j,i,k,1) = (0.5*(1-0.5)*sin(2.0*pi*i*4/N)); %intial cases for laminar, potential fields
        w(j,i,k,2) =-(0.5*(1-0.5)*sin(2.0*pi*i*4/N));
        w(j,i,k,3) = w(j,i,k,1);
        end
    end
end
p  = (w(:,:,:,1)+w(:,:,:,2)+w(:,:,:,3))/3; %p is a N X N matrix

x(1)=0; 
for i=2:N
    x(i) = x(i-1) + dr;
end

%w   = randn(L,L,3);
q   = zeros(N,N,N,NPF);  %IT SHOULD ONLY NE N X N MATRIX
q_i = zeros(N,N,N,NPF);
qg   = zeros(N,N,N,NPG);
qg_i = zeros(N,N,N,NPG);
rho = zeros(N,N,N,3);
rhog= zeros(N,N,N,1);
f   = [ones(1,NSFA+1) 2*ones(1,NSFB)];
f_i = [2*ones(1,NSFB+1) ones(1,NSFA)];
%k^2

%CHANGE THIS TO 3D
L=N;
kk = zeros(L,L,L);
for y=1:L
tk1  = (1:L/2)'*ones(1,L/2);
tk0  = tk1-1;
tk0  = tk0.*tk0;
tk1  = tk1.*tk1;
tk2  = tk0+tk1';
tk3  = tk1+tk1';
tk4=transpose(tk2);
if (y<=(L/2)+1)
   hh=((y-1)*(y-1)).*ones(L,L);
   kk(:,:,y)   =hh(:,:)+ [ tk0+tk0', tk2(:,L/2:-1:1);
    tk4(L/2:-1:1,:),tk3(L/2:-1:1,L/2:-1:1)];
   
else if(y>(L/2)+1)
        hh=((L-y+1)*(L-y+1)).*ones(L,L);
kk(:,:,y)   =hh(:,:)+ [ tk0+tk0', tk2(:,L/2:-1:1);
    tk4(L/2:-1:1,:),tk3(L/2:-1:1,L/2:-1:1)];
 
    end
end    
end




ffe = fopen('fe.dat','wt+');
fprintf(ffe,'%%\n');
fclose(ffe);
        
counter=1;
ft_old=100.0;

for i=1:900
    %pseudo-spectral algorithm
    q(:,:,:,1)=1;   
    q_i(:,:,:,1)=1; 
    for s=2:NPF
        d   = exp(-0.5*ds*w(:,:,:,f(s))).*q(:,:,:,s-1);
        d_i = exp(-0.5*ds*w(:,:,:,f_i(s))).*q_i(:,:,:,s-1);
        fd   = fftn(d);
        fd_i = fftn(d_i);
        tmp  = exp(-4*pi*pi*ds*kk/V);
        fd   = fd.*tmp; %FORMULA FOR FOURIER TRANSFORMATION
        fd_i = fd_i.*tmp;
        q(:,:,:,s)   = exp(-0.5*ds*w(:,:,:,f(s))).*real(ifftn(fd));
        q_i(:,:,:,s) = exp(-0.5*ds*w(:,:,:,f_i(s))).*real(ifftn(fd_i));
    end

    qg(:,:,:,1)=1;
    for s=2:NPG
        dg   = exp(-0.5*dsg*beta*w(:,:,:,1)).*qg(:,:,:,s-1);
        fdg  = fftn(dg);
        tmp  = exp(-4*pi*pi*dsg*beta*kk/V);
        fdg  = fdg.*tmp;
        qg(:,:,:,s)  = exp(-0.5*dsg*beta*w(:,:,:,1)).*real(ifftn(fdg));
    end

    qg_i(:,:,:,1)=exp(-alpha*w(:,:,:,3)).*(qg(:,:,:,NPG).^(sigma-1));
    for s=2:NPG
        dg_i = exp(-0.5*dsg*beta*w(:,:,:,1)).*qg_i(:,:,:,s-1);
        fdg_i= fftn(dg_i);
        tmp  = exp(-4*pi*pi*dsg*beta*kk/V);
        fdg_i= fdg_i.*tmp;
        qg_i(:,:,:,s)  = exp(-0.5*dsg*beta*w(:,:,:,1)).*real(ifftn(fdg_i));
    end
    
    sumphiold = rho(:,:,:,1) + rho(:,:,:,2) + rho(:,:,:,3);
    
    %densify function
    %p = (w(:,:,1)+w(:,:,2)+w(:,:,3))/3;
    %p = p-mean(mean(p)); %shift
    
    for j = 1:N
        for k = 1:N
            for l= 1:N
            a1 = q(j,k,l,:).*q_i(j,k,l,NPF:-1:1);
            rho(j,k,l,1) = simpson1d(a1(1:NSFA+1),ds);
            rho(j,k,l,2) = simpson1d(a1(NSFA+1:NPF),ds);
            a1g= qg(j,k,l,:).*qg_i(j,k,l,NPG:-1:1);
            rhog(j,k,l,1)= simpson1d(a1g(1:NPG),dsg);
            end
        end
    end
    a2  = q(:,:,:,2).*q_i(:,:,:,NPF-1);
    a2A = rho(:,:,:,1);
    a2B = rho(:,:,:,2);
      Qc = simpson3d(a2,dr,dr,dr);
    pfA = simpson3d(a2A,dr,dr,dr);
    pfB = simpson3d(a2B,dr,dr,dr);
    rho(:,:,:,1) = rho(:,:,:,1)*V*phicb*fA/pfA;
    rho(:,:,:,2) = rho(:,:,:,2)*V*phicb*fB/pfB;
    
    a3g = zeros(N,N,N);
    for j = 1:N
        for k = 1:N
            for l=1:N
            a3g(j,k,l) = exp(-alpha*w(j,k,l,3))*(qg(j,k,l,NPG)^sigma);
            end
        end 
    end
%    a3g(:,:) = exp(-alpha*w(:,:,3)).*(qg(:,:,NPG).^sigma);
    Qg = simpson3d(a3g,dr,dr,dr);
    
    a2g = rhog(:,:,:,1);
    pfg = simpson3d(a2g,dr,dr,dr);
    rhog(:,:,:,1)  = rhog(:,:,:,1)*(V*phigpb/pfg)*(sigma*beta/(alpha + sigma*beta));
    
    rho(:,:,:,1) = rho(:,:,:,1) + rhog(:,:,:,1);
    
    rho(:,:,:,3) = exp(-alpha*w(:,:,:,3)).*(qg(:,:,:,NPG).^sigma);
    rho(:,:,:,3) = rho(:,:,:,3)*(V*phigpb/Qg)*(alpha/(alpha + sigma*beta));
    
    sumphinew = rho(:,:,:,1) + rho(:,:,:,2) + rho(:,:,:,3);
    
    %free energy
    tmp = ChiAB*NP*rho(:,:,:,1).*rho(:,:,:,2)+...
          ChiAP*NP*rho(:,:,:,1).*rho(:,:,:,3)+...
          ChiBP*NP*rho(:,:,:,2).*rho(:,:,:,3)-...
          w(:,:,:,1).*rho(:,:,:,1)-...
          w(:,:,:,2).*rho(:,:,:,2)-...
          w(:,:,:,3).*rho(:,:,:,3)-...
          p.*(ones(N,N,N)-rho(:,:,:,1)-rho(:,:,:,2)-rho(:,:,:,3))*0; % MULTIPLIED BY ZERO???
    fe = simpson3d(tmp,dr,dr,dr)/V;
    fs = -phicb*log(Qc/(V*phicb)) - ...
          phigpb*log(Qg/(V*phigpb))/(alpha+sigma*beta);
    ft = fe + fs;
    fprintf(1, 'i = %d, fe=%f, fs=%f, ft=%f\n', i, fe, fs, ft);
    ft_new=ft;
    
    %iterate
    w(:,:,:,1) = w(:,:,:,1) + lambda*(ChiAB*NP*(rho(:,:,:,2)-fB*0)+...
               ChiAP*NP*rho(:,:,:,3) + p - w(:,:,:,1)) -...
               eps*(sumphiold - sumphinew);
    w(:,:,:,2) = w(:,:,:,2) + lambda*(ChiAB*NP*(rho(:,:,:,1)-fA*0)+...
               ChiBP*NP*rho(:,:,:,3) + p - w(:,:,:,2)) -...
               eps*(sumphiold - sumphinew);  
    w(:,:,:,3) = w(:,:,:,3) + lambda*(ChiAP*NP*(rho(:,:,:,1)-fA*0)+...
               ChiBP*NP*(rho(:,:,:,2)-fB*0) + p - w(:,:,:,3)) -...
               eps*(sumphiold - sumphinew);  
    
    p = p - lambda*(1.0 - (rho(:,:,:,1) + rho(:,:,:,2) + rho(:,:,:,3)));
    p = p - mean(mean(mean(p))); %shift
    
    %image((rho(:,:,1)+0.00)*100);
    %colorbar();
    %drawnow
    
     t = 1:L;
scales = 1:L;
x = 1:L;
[T,SCALES,X] = meshgrid(t,scales,x);
tslice = []; 
scalesslice = []; 
xslice = 1:L;
surfHandles = slice(T,SCALES,X,rho(:,:,:,2),tslice,scalesslice,xslice);
set(surfHandles,'FaceAlpha',0.4,'EdgeAlpha',0.1)

   colormap(jet);
    drawnow
    
    if (rem(i,1)==0)
        ffe = fopen('fe.dat','at+');
        fprintf(ffe,'%d\t%.8f\t%.8f\t%.8f\n',i,fe,fs,ft);
        fclose(ffe);
    end 
    
    if (rem(i,10)==0)
        filename = sprintf('phi%.4d.dat', counter);
        fid = fopen(filename,'wt+');
        for j = 1:N
            for k=1:N
                for l=1:N
            fprintf(fid,'%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t %.8f\n',dr*j,dr*k,dr*l,...
                rho(k,j,l,1),rho(k,j,l,2),rho(k,j,l,3));
                end
            end
        end
        fclose(fid);
        counter=counter+1;
    end
    
    ft_diff=abs(ft_new - ft_old);
    if(ft_diff <= 10^-07 || i > 10000)
        break
    else
        ft_old=ft_new;
    end
end