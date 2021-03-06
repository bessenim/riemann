format long;
%% Ѕлок дл€ описани€ начальных данных:
%«асекаем врем€ работы программы:
tic;

sdim = 2;                                                                   %–азмерность по пространственным переменым

%   «адаем параметры задачи
la = 2; mu = 1; ro = 9;                                                     % параметры среды
A1_m_n = ...
[ 0         0   0       -(la+2*mu)  0   ;...
  0         0   0       -la         0   ;...
  0         0   0       0           -mu ;...
  -1/ro     0   0       0           0   ;...
  0         0   -1/ro   0           0   ];

A2_m_n = ...
[ 0 0       0       0   -la         ;...
  0 0       0       0   -(la+2*mu)  ;...
  0 0       0       -mu 0           ;...
  0 0       -1/ro   0   0           ;...
  0 -1/ro   0       0   0           ];

sn = size(A1_m_n,1);    sk1 = sn;   sk2 = sn;

omega = 10;                                                                  %частота действующей силы
% %¬ременные параметры:
T = 1.2;
Ninput= 21;

%   «адаем область, в которой ищем решение
d1min = -0.5; d1max = 0.5;                                                    %Ћева€ и права€ границы области
d2min = -0.5; d2max = 0.5;                                                    %Ћева€ и права€ границы области

%   «адаем пор€док интерпол€ционных полиномов по пространству
K1max = 1;  
K2max = 1;  

%   «адаем сетку по пространственным переменным x_ii1 и x_ii2
II1max = 2*10+1;  I1max = II1max-1; si1 = I1max;                            % оличество узлов сетки-всегда нечетное
x_ii1 = d1min:((d1max-d1min)/si1):d1max;
h_i1 = x_ii1(2:end)-x_ii1(1:(end-1));

II2max = 2*10+1;  I2max = II2max-1; si2 = I2max;                            % оличество узлов сетки-всегда нечетное
x_ii2 = d2min:((d2max-d2min)/si2):d2max;
h_i2 = x_ii2(2:end)-x_ii2(1:(end-1));

%%   ¬спомагательные процедуры
%   l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 
sl1 = K1max + 1;
ksi_l1 = zeros(1, sl1);
if sl1==1
   ksi_l1 = 0;
else
    ksi_l1 = -1:(2/(sl1-1)):1;
%     ksi_l1(1, :) = -cos(pi/(2*(K1max+1))+pi*(0:1:K1max)/(K1max+1));
end
x_i1_l1 = zeros( si1, sl1 );
for i1 = 1:si1
    x_i1_l1(i1,:) = 0.5*(ksi_l1*(x_ii1(i1+1)-x_ii1(i1))+(x_ii1(i1+1)+x_ii1(i1)));
end

% ksi_l1(1, :) = -cos(pi/(2*(K1max+1))+pi*(0:1:K1max)/(K1max+1));
ksi_l1_k1 = zeros(sl1,sl1);
for k1 = 0:K1max
    ksi_l1_k1(:, k1+1) = ksi_l1(1, :).^(K1max-k1);
end

pbase_l1_k1 = zeros(sl1, sl1);
pbase_l1_k1 = (ksi_l1_k1\eye(sl1))';                                          %   pbase_l1_k1 = PolyInterpChebNode_k1_l1

%       —троим вспомагательные массивы
sp1=si1*(sl1-1)+1;
V_l1_i1_p1=zeros(sl1,si1,sp1);
invV_p1_l1_i1=zeros(sp1,sl1,si1);
for i1=1:si1
    for l1=1:(sl1-1)
        p1=(i1-1)*(sl1-1)+l1;
        V_l1_i1_p1(l1,i1,p1)=1;
    end
    l1=sl1;
    p1=(i1)*(sl1-1)+1;
    V_l1_i1_p1(l1,i1,p1)=1;
end

for i1=1:si1
    for l1=1:(sl1-1)
        p1=(i1-1)*(sl1-1)+l1;
        invV_p1_l1_i1(p1,l1,i1)=1;
    end
end
invV_p1_l1_i1(sp1,sl1,si1)=1;

V_i1_l1_p1 = permute(V_l1_i1_p1, [ 2 1 3 ]);
invV_p1_i1_l1 = permute(invV_p1_l1_i1, [ 1 3 2 ]);
clear V_l1_i1_p1 invV_p1_l1_i1;

x_p1 = reshape(invV_p1_i1_l1, [ sp1 si1*sl1 ])*reshape(x_i1_l1, [ si1*sl1 1 ]);

sl2 = K2max + 1;
ksi_l2 = zeros(1, sl2);
if sl2==1
   ksi_l2 = 0;
else
    ksi_l2 = -1:(2/(sl2-1)):1;
%     ksi_l2(1, :)=-cos(pi/(2*(K2max+1))+pi*(0:1:K2max)/(K2max+1));
end
x_i2_l2 = zeros( si2, sl2 );
for i2 = 1:si2
    x_i2_l2(i2,:) = 0.5*(ksi_l2*(x_ii2(i2+1)-x_ii2(i2))+(x_ii2(i2+1)+x_ii2(i2)));
end

% ksi_l2(1, :)=-cos(pi/(2*(K2max+1))+pi*(0:1:K2max)/(K2max+1));
ksi_l2_k2 = zeros(sl2,sl2);
for k2 = 0:K2max
    ksi_l2_k2(:, k2+1) = ksi_l2(1, :).^(K2max-k2);
end

pbase_l2_k2 = zeros(sl2, sl2);
pbase_l2_k2 = (ksi_l2_k2\eye(sl2))';                                          %   pbase_l2_k2 = PolyInterpChebNode_k2_l2

%       —троим вспомагательные массивы
sp2=si2*(sl2-1)+1;
V_l2_i2_p2=zeros(sl2,si2,sp2);
invV_p2_l2_i2=zeros(sp2,sl2,si2);
for i2=1:si2
    for l2=1:(sl2-1)
        p2=(i2-1)*(sl2-1)+l2;
        V_l2_i2_p2(l2,i2,p2)=1;
    end
    l2=sl2;
    p2=(i2)*(sl2-1)+1;
    V_l2_i2_p2(l2,i2,p2)=1;
end

for i2=1:si2
    for l2=1:(sl2-1)
        p2=(i2-1)*(sl2-1)+l2;
        invV_p2_l2_i2(p2,l2,i2)=1;
    end
end
invV_p2_l2_i2(sp2,sl2,si2)=1;

V_i2_l2_p2 = permute(V_l2_i2_p2, [ 2 1 3 ]);
invV_p2_i2_l2 = permute(invV_p2_l2_i2, [ 1 3 2 ]);
clear V_l2_i2_p2 invV_p2_l2_i2;

x_p2 = reshape(invV_p2_i2_l2, [ sp2 si2*sl2 ])*reshape(x_i2_l2, [ si2*sl2 1 ]);

x_i_p1_p2 = zeros(sdim,sp1,sp2);
for p1 = 1:sp1
    x_i_p1_p2(1,p1,:) = x_p1(p1);
end
for p2 = 1:sp2
    x_i_p1_p2(2,:,p2) = x_p2(p2);
end
x_p1_p2_i = permute(x_i_p1_p2, [ 2 3 1 ]);

% ћассивы C_k1_m_n, C_k2_m_n
[R1_m_n,LAMBD1_m_n] = eig(A1_m_n);
lambd1 = diag(LAMBD1_m_n);  
% indlambd1lt0 = find(lambd1<0);  %   lengthindlambd1lt0 = length(indlambd1lt0);
indlambd1le0 = find(lambd1<=0); %   lengthindlambd1le0 = length(indlambd1le0);
% indlambd1gt0 = find(lambd1>0);  %   lengthindlambd1gt0 = length(indlambd1gt0);
indlambd1ge0 = find(lambd1>=0); %   lengthindlambd1ge0 = length(indlambd1ge0);

L1_m_n = inv(R1_m_n);
D_k1_m_n = zeros(sn,sn,sn);
for k1=1:sn
    D_k1_m_n(k1,k1,k1) = 1;
end
C_k1_m_n = zeros(sn,sn,sn);
for k1=1:sn
    C_k1_m_n(k1,:,:) = R1_m_n*reshape(D_k1_m_n(k1,:,:), [ sn sn ])*L1_m_n;
end
clear A1_m_n R1_m_n LAMBD1_m_n L1_m_n D_k1_m_n ;

[R2_m_n,LAMBD2_m_n] = eig(A2_m_n);
lambd2 = diag(LAMBD2_m_n);
% indlambd2lt0 = find(lambd2<0);  
indlambd2le0 = find(lambd2<=0); 
% indlambd2gt0 = find(lambd2>0);  
indlambd2ge0 = find(lambd2>=0); 

L2_m_n = inv(R2_m_n);
L2GE0_m_n = L2_m_n(find(lambd2>=0),:);
D_k2_m_n = zeros(sn,sn,sn);
for k2=1:sn
    D_k2_m_n(k2,k2,k2) = 1;
end
C_k2_m_n = zeros(sn,sn,sn);
for k2=1:sn
    C_k2_m_n(k2,:,:) = R2_m_n*reshape(D_k2_m_n(k2,:,:), [ sn sn ])*L2_m_n;
end
clear A2_m_n R2_m_n LAMBD2_m_n L2_m_n D_k2_m_n ;

deltat = min(min(h_i1)/max(abs(lambd1)), min(h_i2)/max(abs(lambd2))) * 0.25;
Ntau = fix(T/deltat);
deltaNinput = fix((Ntau+1)/(Ninput-1));
%   it = 1:1:(deltaNinput*Ninput);

C_k1_k2_m_n = zeros(sn,sn,sn,sn);
for k1=1:sn
    for k2=1:sn
        C_k1_k2_m_n(k1,k2,:,:) = 0.5* ...
            (  reshape(C_k1_m_n(k1,:,:), [ sn sn ])*reshape(C_k2_m_n(k2,:,:), [ sn sn ]) ...
             + reshape(C_k2_m_n(k2,:,:), [ sn sn ])*reshape(C_k1_m_n(k1,:,:), [ sn sn ]));
    end
end
clear C_k1_m_n C_k2_m_n;
C_k2k1n_m = reshape(permute(C_k1_k2_m_n, [ 2 1 4 3 ]), [ sk2*sk1*sn sn ]);

[numlambd2GE0,ind] = size(L2GE0_m_n);
L2GE0C_k1_k2_m_n = zeros(sn,sn,numlambd2GE0,sn);
for k1 = 1:sn
    for k2 = 1:sn
        L2GE0C_k1_k2_m_n(k1,k2,:,:) = reshape(...
            L2GE0_m_n*reshape(C_k1_k2_m_n(k1,k2,:,:), [ sn sn ]), [ 1 1 numlambd2GE0 sn ]);
    end
end

H_p1w_k1_p1 = zeros(sp1,sn,sp1);
for p1w = 1
    for k1 = indlambd1le0
H_p1w_k1_p1(p1w,k1,p1w)   = polyval(pbase_l1_k1(1,:), (-1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
H_p1w_k1_p1(p1w,k1,p1w+1) = polyval(pbase_l1_k1(2,:), (-1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
    end
end

for p1w = 2:(sp1-1)
    for k1 = indlambd1le0
H_p1w_k1_p1(p1w,k1,p1w)   = polyval(pbase_l1_k1(1,:), (-1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
H_p1w_k1_p1(p1w,k1,p1w+1) = polyval(pbase_l1_k1(2,:), (-1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
    end
    for k1 = indlambd1ge0
H_p1w_k1_p1(p1w,k1,p1w)   = polyval(pbase_l1_k1(2,:), (1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
H_p1w_k1_p1(p1w,k1,p1w-1) = polyval(pbase_l1_k1(1,:), (1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
    end
end

for p1w = sp1
    for k1 = indlambd1ge0
H_p1w_k1_p1(p1w,k1,p1w)   = polyval(pbase_l1_k1(2,:), (1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
H_p1w_k1_p1(p1w,k1,p1w-1) = polyval(pbase_l1_k1(1,:), (1-(2/(h_i1(i1)))*lambd1(k1)*deltat));
    end
end

H_p2w_k2_p2 = zeros(sp2,sn,sp2);
for p2w = 1
    for k2 = indlambd2le0
H_p2w_k2_p2(p2w,k2,p2w)   = polyval(pbase_l2_k2(1,:), (-1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
H_p2w_k2_p2(p2w,k2,p2w+1) = polyval(pbase_l2_k2(2,:), (-1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
    end
end

for p2w = 2:(sp2-1)
    for k2 = indlambd2le0
H_p2w_k2_p2(p2w,k2,p2w)   = polyval(pbase_l2_k2(1,:), (-1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
H_p2w_k2_p2(p2w,k2,p2w+1) = polyval(pbase_l2_k2(2,:), (-1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
    end
    for k2 = indlambd2ge0
H_p2w_k2_p2(p2w,k2,p2w)   = polyval(pbase_l2_k2(2,:), (1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
H_p2w_k2_p2(p2w,k2,p2w-1) = polyval(pbase_l2_k2(1,:), (1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
    end
end

for p2w = sp2
    for k2 = indlambd2ge0
H_p2w_k2_p2(p2w,k2,p2w)   = polyval(pbase_l2_k2(2,:), (1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
H_p2w_k2_p2(p2w,k2,p2w-1) = polyval(pbase_l2_k2(1,:), (1-(2/(h_i2(i2)))*lambd2(k2)*deltat));
    end
end

A1_m_n = zeros(2,sn);
A1_m_n(1,3) = 1;
A1_m_n(2,2) = 1;

A_m_n = cat(1,L2GE0_m_n,A1_m_n);
invA_m_n = inv(A_m_n);

%%  Ќачальные данные u0_p1_p2_n
u0_p1_p2_n = zeros( sp1,sp2,sn );
amplitudeu0 = 0.1;
xstar_i = [0 0];
diam = 0.05;
centerline = 0.1;
radius_p1_p2_i = zeros(sp1,sp2,sdim);
normradius_p1_p2 = zeros(sp1,sp2);
for p1 = 1:sp1
    for p2 = 1:sp2
        radius_p1_p2_i(p1,p2,:) = x_p1_p2_i(p1,p2,:) - reshape(xstar_i, [ 1 1 sdim ]);
        normradius_p1_p2(p1,p2) = norm(reshape(radius_p1_p2_i(p1,p2,:), [ sdim 1 ]),2);
    end
end

[indp1,indp2] = find(abs(normradius_p1_p2-centerline) <= diam);
lengthind = length(indp1);
for point = 1:lengthind
    u0_p1_p2_n(indp1(point),indp2(point),4:5) = amplitudeu0*...
       (radius_p1_p2_i(indp1(point),indp2(point),:)/normradius_p1_p2(indp1(point),indp2(point)))...
       *(1 + cos((pi / diam)*(normradius_p1_p2(indp1(point),indp2(point))-centerline)));
end
% 
% for p2 = 1:((sp2-1)/2 - 1)
%     for p1 = 1:sp1
%         u0_p1_p2_n(p1,p2,1) = 0;
%         u0_p1_p2_n(p1,p2,2) = 0;
%         u0_p1_p2_n(p1,p2,3) = 0;
%         u0_p1_p2_n(p1,p2,4) = 0;
%         u0_p1_p2_n(p1,p2,5) = 0;
%     end
% end
% for p2 = ((sp2-1)/2)
%     for p1 = 1:sp1
%         u0_p1_p2_n(p1,p2,1) = 0;
%         u0_p1_p2_n(p1,p2,2) = 0;
%         u0_p1_p2_n(p1,p2,3) = 0.948683298050514/2;
%         u0_p1_p2_n(p1,p2,4) = 0;
%         u0_p1_p2_n(p1,p2,5) = -0.316227766016838/2;
%     end
% end
% for p2 = ((sp2-1)/2 + 1):sp2
%     for p1 = 1:sp1
%         u0_p1_p2_n(p1,p2,1) = 0;
%         u0_p1_p2_n(p1,p2,2) = 0;
%         u0_p1_p2_n(p1,p2,3) = 0.948683298050514;
%         u0_p1_p2_n(p1,p2,4) = 0;
%         u0_p1_p2_n(p1,p2,5) = -0.316227766016838;
%     end
% end

%%   ƒелаем шаги по времени
sjt = Ninput;
T_jt = 0:deltat*deltaNinput:(sjt-1)*deltat*deltaNinput;
% [per, sjt] = size(T_jt);   clear per;

u_p1_p2_n_jt = zeros(sp1,sp2,sn,Ninput);
uprev_p1_p2_n = zeros(sp1,sp2,sn);
unext_p1_p2_n = zeros(sp1,sp2,sn);
u_p1_p2_n_jt(:,:,:,1) = u0_p1_p2_n;
uprev_p1_p2_n(:,:,:) = u_p1_p2_n_jt(:,:,:,1);
indk1 = 1:sk1;  indk2 = 1:sk2;


amplitude_p1w_sp2 = zeros(sp1,1);   %% ИСТОЧНИКИ
%amplitude_p1w_sp2((sp1+1)/2-1, 1) = 10;    
amplitude_p1w_sp2((sp1+1)/2, 1) = 0;
%amplitude_p1w_sp2((sp1+1)/2+1, 1) = 10;
%amplitude_p1w_sp2(round((sp1+1)/3)-1, 1) = 10;    
%amplitude_p1w_sp2(round((sp1+1)/3), 1) = 10;
%amplitude_p1w_sp2(round((sp1+1)/3)+1, 1) = 10;
%amplitude_p1w_sp2(round(2*(sp1+1)/3)-1, 1) = 10;    
%amplitude_p1w_sp2(round(2*(sp1+1)/3), 1) = 10;
%amplitude_p1w_sp2(round(2*(sp1+1)/3)+1, 1) = 10;
        
ind_p1w_p1 = zeros(sp1,3);  lengthind_p1w = zeros(sp1,1);
p1w = 1;    ind_p1w_p1(p1w,1) = p1w;    ind_p1w_p1(p1w,2) = p1w+1;  lengthind_p1w(p1w) = 2;
for p1w = 2:(sp1-1)
    ind_p1w_p1(p1w,1) = p1w-1;  ind_p1w_p1(p1w,2) = p1w;  ind_p1w_p1(p1w,3) = p1w+1;  lengthind_p1w(p1w) = 3;
end
p1w = sp1;    ind_p1w_p1(p1w,1) = p1w-1;    ind_p1w_p1(p1w,2) = p1w;  lengthind_p1w(p1w) = 2;

ind_p2w_p2 = zeros(sp2,3);  lengthind_p2w = zeros(sp2,1);
p2w = 1;    ind_p2w_p2(p2w,1) = p2w;    ind_p2w_p2(p2w,2) = p2w+1;  lengthind_p2w(p2w) = 2;
for p2w = 2:(sp2-1)
    ind_p2w_p2(p2w,1) = p2w-1;  ind_p2w_p2(p2w,2) = p2w;  ind_p2w_p2(p2w,3) = p2w+1;  lengthind_p2w(p2w) = 3;
end
p2w = sp2;    ind_p2w_p2(p2w,1) = p2w-1;    ind_p2w_p2(p2w,2) = p2w;  lengthind_p2w(p2w) = 2;


C_k2_k1_n_m = permute(C_k1_k2_m_n, [ 2 1 4 3 ]);
C_k2k1n_m = reshape(C_k2_k1_n_m, [ sn*sn*sn sn ]);
for jt=2:sjt
    for point = 2:deltaNinput+1
        for p1w = 1:sp1
W_k1_p2_n = reshape(reshape(H_p1w_k1_p1(p1w,:,ind_p1w_p1(p1w,1:lengthind_p1w(p1w))), [ sk1 lengthind_p1w(p1w) ])...
                   *reshape(uprev_p1_p2_n(ind_p1w_p1(p1w,1:lengthind_p1w(p1w)),:,:), [ lengthind_p1w(p1w) sp2*sn ]), [ sk1 sp2 sn ]);
W_p2_k1_n = permute(W_k1_p2_n, [ 2 1 3 ]);
             for p2w = 1:(sp2-1)
W_k2k1n = reshape(reshape(H_p2w_k2_p2(p2w,:,ind_p2w_p2(p2w,1:lengthind_p2w(p2w))), [ sk2 lengthind_p2w(p2w) ])...
                   *reshape(W_p2_k1_n(ind_p2w_p2(p2w,1:lengthind_p2w(p2w)),:,:), [ lengthind_p2w(p2w) sk1*sn ]), [ 1 sk2*sk1*sn ]);
unext_p1_p2_n(p1w,p2w,:) = reshape(W_k2k1n*C_k2k1n_m, [ 1 1 sn ]);               
%                 unext_p1_p2_n(p1w,p2w,:)
             end
            for p2w = sp2
W_k2_k1_n = reshape(reshape(H_p2w_k2_p2(p2w,:,ind_p2w_p2(p2w,1:lengthind_p2w(p2w))), [ sk2 lengthind_p2w(p2w) ])...
                   *reshape(W_p2_k1_n(ind_p2w_p2(p2w,1:lengthind_p2w(p2w)),:,:), [ lengthind_p2w(p2w) sk1*sn ]), [ sk2 sk1 sn ]);
                W_m = zeros(sn,1);
                for k1 = indk1
                    for k2 = indk2(lambd2>=0)
                        W_m(:,1) = W_m(:,1) + ...
                            reshape(C_k1_k2_m_n(k1,k2,:,:), [ sn sn ])...
                           *reshape(W_k2_k1_n(k2,k1,:), [ sn 1 ]);
                    end
                end
                r1_m = L2GE0_m_n*W_m;
                t = (jt-1)*deltat*deltaNinput + point*deltat;
                r2_m(1,1) = 0; r2_m(2,1) = -amplitude_p1w_sp2(p1w,1)*sin(omega*t);
                r_m = cat(1,r1_m,r2_m); 
                unext_p1_p2_n(p1w,p2w,:) = invA_m_n*r_m;
            end
        end
        uprev_p1_p2_n = unext_p1_p2_n;
    end
    u_p1_p2_n_jt(:,:,:,jt) = unext_p1_p2_n(:,:,:);
end

K_p1_p2_jt = zeros(sp1,sp2,sjt);
for p1 = 1:sp1
    for p2 = 1:sp2
        for jt = 1:sjt
            %K_p1_p2_jt(p1,p2,jt) = u_p1_p2_n_jt(p1,p2,4,jt);
            K_p1_p2_jt(p1,p2,jt) = sqrt((u_p1_p2_n_jt(p1,p2,4,jt).^2 + u_p1_p2_n_jt(p1,p2,5,jt).^2))*1;
        end
    end
end

toc;

%%  √рафический вывод результатов
x1min = min(x_p1);                       x1max = max(x_p1);
x2min = min(x_p2);                       x2max = max(x_p2);

Kmin = min(min(min(K_p1_p2_jt, [], 1), [], 2), [], 3);  
Kmax = max(max(max(K_p1_p2_jt, [], 1), [], 2), [], 3);

% [X1,X2] = meshgrid(x_p1,x_p2);

[per, sit] = size(T_jt);    clear per;
[Y, X] = meshgrid(x_p2, x_p1);
        for it=1:1:sit
            %surface = reshape(K_p1_p2_jt(:,:,it), [ sp1 sp2 ]);
            %vtkwrite('surf.vtk','structured_grid',x_p1,x_p2,'scalars', surface)
            %surf(x_p2,x_p1,reshape(K_p1_p2_jt(:,:,it), [ sp1 sp2 ]));
            %axis([ x2min x2max x1min x1max Kmin Kmax]);
            %pause(1);
            %F(it) = getframe; %#ok<AGROW>
            %tic
            name = sprintf('res%d.vtk', it);
            f = fopen(name, 'wb');
            fprintf(f, '# vtk DataFile Version 3.0\n');
            fprintf(f, 'Exported from MATLAB\n'); % Comment string
            fprintf(f, 'BINARY\n');
            fprintf(f, 'DATASET STRUCTURED_GRID\n');
            fprintf(f, 'DIMENSIONS %d %d 1\n', sp1, sp2);
            fprintf(f, 'POINTS %d float\n', sp1 * sp2);
            R = zeros(3, sp1, sp2);
            R(1, :, :) = Y;
            R(2, :, :) = X;
            R(3, :, :) = K_p1_p2_jt(:,:,it);
            w = typecast(swapbytes(single(R(:))), 'uint8');
            fwrite(f, w);
            fprintf(f, 'CELL_DATA %d\n', (sp1-1) * (sp2-1));
            % No cell data
            fprintf(f, 'POINT_DATA %d\n', sp1 * sp2);
            fprintf(f, 'SCALARS z float\nLOOKUP_TABLE default\n');
            w = typecast(swapbytes(single(reshape(K_p1_p2_jt(:,:,it),1, []))), 'uint8');
            fwrite(f, w);
            fclose(f);
            %toc
        end
        
        %movie(F, 3, 5);

%end
K2_sq = zeros(11);
for i = 1:1:11
    for j = 1:1:11
        K2_sq(i, j) = K_p1_p2_jt(2*i-1,2*j-1,10);
    end
end

