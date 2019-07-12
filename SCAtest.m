% user_K=4,userANt_Nr=4,tansAnt_Nt=8,timelength_T=10,QAM_orderM=16
...ser_Sigma=0.01(for all),delata=1(for all)
% transignal_S:
K=4;
Nr=4;
Nt=8;
T=10;
M=16;
rng(100);
S=qammod(randi(16,K,T)-1,16);
S1=[real(S);imag(S)];
% alpha and c
alpha=zeros(2*K,T);
c=zeros(2*K,T);
alpha(:,:)=1/sqrt(2)*qfuncinv(0.01/2);
c(:,:)=1/sqrt(2)*qfuncinv(0.01/2);
Index_pM=find(S1==sqrt(M)-1);
Index_mM=find(S1==-sqrt(M)+1);
alpha(Index_pM)=1/sqrt(2)*qfuncinv(0.01);
alpha(Index_mM)=-inf;
c(Index_pM)=-inf;
c(Index_mM)=1/sqrt(2)*qfuncinv(0.01);
% channel
Ht=zeros(Nr,Nt,K);
for i=1:K
    Ht(:,:,i)=randn(Nr,Nt)+1i*randn(Nr,Nt)
end
%one time cvx
cvx_begin quiet
    variable d(K,1) nonnegative;
    variable u(K,T) complex;
    variable w(Nr,K) complex;
    expression f(T,1);
    expression G(K,Nt);
    for i=1:K
        G(i,:)=w(:,i)'*Ht(:,:,i);
    end
    for j=1:T
        f(j,1)=(diag(d)*S(:,j)+u(:,j))'*inv(G*G')*(diag(d)*S(:,j)+u(:,j));
    end
    minimize sum(f)
        subject to
           for k=1:T
               -d+alpha(1:K,k)<=real(u);
               real(u)<=d-c(1:K,k);
               -d+alpha(K+1:end,k)<=imag(u);
               imag(u)<=d-c(K+1:end,k);
           end
           for l=1:K
               norm(w(:,l))^2<=1;
           end
 cvx_end
               
               
        
