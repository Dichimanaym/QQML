% 3D double linear iteration; hard wall or periodic boundary condiction
clear;close all
tic
p=0; %periodic: p=1, hardwall: p=0;
Lit=0;
for L=6:4:30,
    Lit=Lit+1
    Ln(Lit)=L;
    vait=0;
    col=jet(11);
    clear Gn1
    for va=15:1:20,
        vait=vait+1
        van(vait)=va;
X=L;Y=L;Z=L; % dimensions
%va=15.0; % disorder strength
% hopping elements of the slices
HX=diag(ones(X-1,1),1)+diag(ones(X-1,1),-1);HX(1,X)=p;HX(X,1)=p;
HY=diag(ones(Y-1,1),1)+diag(ones(Y-1,1),-1);HY(1,Y)=p;HY(Y,1)=p;
H1=kron(HX,eye(Y))+kron(eye(X),HY);
H1=gpuArray(H1);
En=[-.1:.005:.1];
for Eit=1:length(En) % energy iteration
E=En(Eit);
XL=X; % width of the lead
YL=Y; % thickness of the lead
T=eye(X*Y);% connecting two slices (T_{n,n+1})
T=gpuArray(T);
    %% Random lead
HX=diag(ones(XL-1,1),1)+diag(ones(XL-1,1),-1);HX(1,XL)=p;HX(XL,1)=p;
HY=diag(ones(YL-1,1),1)+diag(ones(YL-1,1),-1);HY(1,YL)=p;HY(YL,1)=p;
VL=(rand(XL*YL,1)-.5)*va;
HL=kron(HX,eye(YL))+kron(eye(XL),HY)+diag(VL);% Hamiltonian of the lead
%HL=gpuArray(HL);
[V2,E2]=eig(HL);E2=diag(E2);
k1=real(acos((E-E2)/2))+i*abs(imag(acos((E-E2)/2)));
for x1=1:X*Y, for x2=1:X*Y, Self0R(x1,x2)=sum(V2(x1,:).*conj(V2(x2,:)).*transpose(exp(i*k1)));      
end;end;
Self0R=gpuArray(Self0R);
Self0L=Self0R;
for avv=1:10, % disorder average
%%
Gn=Self0R;
U=(rand(X*Y,Z)-.5)*va; % disorder
U=gpuArray(U);
% first iteration from end to start
Emat=(E-1e-10*i)*eye(X*Y);
Emat=gpuArray(Emat);
clear psi LDOS
for n=Z:-1:1, % iteration for nth slice
    Vn=U(:,n);
Hn=diag(Vn)+H1;
Hn=gpuArray(Hn);
Gn=inv(Emat-Hn-T*Gn*T');Gn1(:,:,n)=Gn;
end;
% second iteration from start to end
Gn2=Self0L;
G1n2=T';
psi0=-i*ones(X*Y,1);
for n=1:Z-1, % iteration for nth slice
Vn=U(:,n);
Hn=diag(Vn)+H1;
Gn=inv(Emat-Hn-T*Gn1(:,:,n+1)*T'-T'*Gn2*T);
Gn2=inv(Emat-Hn-T'*Gn2*T);
G1n=G1n2*T*Gn;
G1n2=G1n2*T*Gn2;
psi(:,n)=-G1n'*psi0;
LDOS(:,n)=imag(diag(Gn));
end;
n=Z; % last slice
Vn=U(:,n);
Hn=diag(Vn)+H1;
Gn=inv(Emat-Hn-Self0R-Gn2);
Gn2=inv(Emat-Hn-Gn2);
G1n=G1n2*Gn;
G1n2=G1n2*Gn2;
psi(:,n)=-G1n'*psi0; % The full solution
LDOS(:,n)=imag(diag(Gn));
GammaR=i*(Self0R-Self0R');
GammaL=i*(Self0L-Self0L');
Trans(Eit,avv,vait,Lit)=trace((GammaR*G1n)*(GammaL*G1n'));
end;
toc
end;
figure(1)
semilogy(En,abs(mean(Trans(:,:,vait,Lit),2)),'color',col(vait,:));hold on
axis([min(En) max(En) .01 10])
    end;
save L6to30va15to20av100En0b En Trans Ln van
end;
% plot(En,abs(Trans))
% plot(En,exp(mean(log(abs(Trans)),2)))
% plot(En,abs(mean(Trans,2)))

