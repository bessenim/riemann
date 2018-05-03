% Ётот файл получен из global_init2_1_3
% јдаптирован дл€ решени€ системы уравнений динамической упругости
% многомерной модификацией сеточно-характеристического метода, 
% использующей фундаментальное решение оператора задачи.

% (c) —калько ёрий, ћ‘“» 2017

%function grid_charact()
format long;
%% Ѕлок дл€ описани€ начальных данных:
%«асекаем врем€ работы программы:
tic;

sdim = 2;                                                                   %–азмерность по пространственным переменым

%   «адаем параметры задачи
la = 2000000; mu = 1000000; ro = 9;                                                     % параметры среды
sn = 5;
A_i_m_n = zeros(sdim,sn,sn);
A_i_m_n(1,:,:) = ...
[ 0         0   0       -(la+2*mu)  0   ;...
  0         0   0       -la         0   ;...
  0         0   0       0           -mu ;...
  -1/ro     0   0       0           0   ;...
  0         0   -1/ro   0           0   ];

A_i_m_n(2,:,:) = ...
[ 0 0       0       0   -la         ;...
  0 0       0       0   -(la+2*mu)  ;...
  0 0       0       -mu 0           ;...
  0 0       -1/ro   0   0           ;...
  0 -1/ro   0       0   0           ];

sg = 9;                                                                     % количество регионов

omega = 10;                                                                  %частота действующей силы
% %¬ременные параметры:
T = 1.1;
Ninput= 21;

%   «адаем область, в которой ищем решение
dmin_i = zeros(1,sdim); dmax_i = zeros(1,sdim);
dmin_i(1) = -30;  dmax_i(1) = 30;
dmin_i(2) = -600;  dmax_i(2) = 0;

dmin_g_i = zeros(sg,sdim); dmax_g_i = zeros(sg,sdim);
dmin_g_i(1,1) = dmin_i(1); dmax_g_i(1,1) = dmin_i(1)/2;
    dmin_g_i(2,1) = dmin_i(1)/2; dmax_g_i(2,1) = dmax_i(1)/2;
        dmin_g_i(3,1) = dmax_i(1)/2; dmax_g_i(3,1) = dmax_i(1);            % Ћева€ и права€ границы области
            dmin_g_i(4,1) = dmin_i(1); dmax_g_i(4,1) = dmin_i(1)/2;
                dmin_g_i(5,1) = dmin_i(1)/2; dmax_g_i(5,1) = dmax_i(1)/2;
                    dmin_g_i(6,1) = dmax_i(1)/2; dmax_g_i(6,1) = dmax_i(1);  
                        dmin_g_i(7,1) = dmin_i(1); dmax_g_i(7,1) = dmin_i(1)/2;
                            dmin_g_i(8,1) = dmin_i(1)/2; dmax_g_i(8,1) = dmax_i(1)/2;
                                   dmin_g_i(9,1) = dmax_i(1)/2; dmax_g_i(9,1) = dmax_i(1);  
dmin_g_i(1,2) = dmin_i(2)/3; dmax_g_i(1,2) = dmax_i(2);
    dmin_g_i(2,2) = dmin_i(2)/3; dmax_g_i(2,2) = dmax_i(2);
        dmin_g_i(3,2) = dmin_i(2)/3; dmax_g_i(3,2) = dmax_i(2);            % верхние и нижние границы области
            dmin_g_i(4,2) = 2*dmin_i(2)/3; dmax_g_i(4,2) = dmin_i(2)/3;
                dmin_g_i(5,2) = 2*dmin_i(2)/3; dmax_g_i(5,2) = dmin_i(2)/3;
                    dmin_g_i(6,2) = 2*dmin_i(2)/3; dmax_g_i(6,2) = dmin_i(2)/3;  
                        dmin_g_i(7,2) = dmin_i(2); dmax_g_i(7,2) = 2*dmin_i(2)/3;
                            dmin_g_i(8,2) = dmin_i(2); dmax_g_i(8,2) = 2*dmin_i(2)/3;
                                   dmin_g_i(9,2) = dmin_i(2); dmax_g_i(9,2) = 2*dmin_i(2)/3;  
    
%   «адаем пор€док интерпол€ционных полиномов по пространству
Kmax = 1;  

%   «адаем сетку по пространственным переменным x_ii1 и x_ii2
sp_i = zeros(1,sdim);   sp_i(1) = 2*30+1;  sp_i(2) = 2*300+1;              %  оличество узлов сетки-всегда нечетное
sp_g_i = zeros(sg,sdim);
sp_g_i(1,1) = (sp_i(1)-1)/4 + 1;
    sp_g_i(2,1) = 2*(sp_i(1)-1)/4 + 1;
        sp_g_i(3,1) = (sp_i(1)-1)/4 + 1;
            sp_g_i(4,1) = (sp_i(1)-1)/4 + 1;
                sp_g_i(5,1) = 2*(sp_i(1)-1)/4 + 1;
                    sp_g_i(6,1) = (sp_i(1)-1)/4 + 1;
                        sp_g_i(7,1) = (sp_i(1)-1)/4 + 1;
                            sp_g_i(8,1) = 2*(sp_i(1)-1)/4 + 1;
                                sp_g_i(9,1) = (sp_i(1)-1)/4 + 1;
sp_g_i(1,2) = (sp_i(2)-1)/3 + 1;
    sp_g_i(2,2) = (sp_i(2)-1)/3 + 1;
        sp_g_i(3,2) = (sp_i(2)-1)/3 + 1;
             sp_g_i(4,2) = (sp_i(2)-1)/3 + 1;
                sp_g_i(5,2) = (sp_i(2)-1)/3 + 1;
                    sp_g_i(6,2) = (sp_i(2)-1)/3 + 1;
                         sp_g_i(7,2) = (sp_i(2)-1)/3 + 1;
                            sp_g_i(8,2) = (sp_i(2)-1)/3 + 1;
                                sp_g_i(9,2) = (sp_i(2)-1)/3 + 1;

for g = 1:sg
    for i = 1:sdim
        x_p{g,i} = dmin_g_i(g,i):((dmax_g_i(g,i)-dmin_g_i(g,i))/(sp_g_i(g,i)-1)):dmax_g_i(g,i);
    end
end
for g = 1:sg
    for i = 1:sdim
        h_g_i_p{g,i} = x_p{g,i}(2:end) - x_p{g,i}(1:(end-1));
    end
end

% x_i_p1_p2 = zeros(sdim,sp1,sp2);
% for p1 = 1:sp1
%     x_i_p1_p2(1,p1,:) = x_p1(p1);
% end
% for p2 = 1:sp2
%     x_i_p1_p2(2,:,p2) = x_p2(p2);
% end
% x_p1_p2_i = permute(x_i_p1_p2, [ 2 3 1 ]);


    
%%   ¬спомагательные процедуры
%   l l l l l l l l l l l l l l l l l l l l l l l 
sl = Kmax + 1;
ksi_l = zeros(1, sl);
if sl==1
   ksi_l = 0;
else
    ksi_l = -1:(2/(sl-1)):1;
%     ksi_l(1, :) = -cos(pi/(2*(Kmax+1))+pi*(0:1:Kmax)/(Kmax+1));
end
ksi_l_k = zeros(sl,sl);
for k = 0:Kmax
    ksi_l_k(:, k+1) = ksi_l(1, :).^(Kmax-k);
end

pbase_l_k = zeros(sl, sl);
pbase_l_k = (ksi_l_k\eye(sl))';                                          %   pbase_l_k = PolyInterpChebNode_k_l

% ћассивы C_i_k_m_n
D_k_m_n = zeros(sn,sn,sn);
for k=1:sn
    D_k_m_n(k,k,k) = 1;
end
C_i_k_m_n = zeros(sdim,sn,sn,sn);

for i = 1:sdim
    [R_i_m_n(i,:,:),LAMBD_i_m_n(i,:,:)] = eig(squeeze(A_i_m_n(i,:,:)));
    lambd_i_k(i,:) = diag(squeeze(LAMBD_i_m_n(i,:,:)));
% %     indlambdlt0_i(i,:) = find(squeeze(lambd_i_k(i,:)<0));   lengthindlambdlt0_i(i) = length(indlambdlt0_i(i,:));
%     indlambdle0_i(i,:) = find(squeeze(lambd_i_k(i,:)<=0));  lengthindlambdle0_i(i) = length(indlambdle0_i(i,:));
% %     indlambdgt0_i(i,:) = find(squeeze(lambd_i_k(i,:)>0));   lengthindlambdgt0_i(i) = length(indlambdgt0_i(i,:));
%     indlambdge0_i(i,:) = find(squeeze(lambd_i_k(i,:)>=0));  lengthindlambdge0_i(i) = length(indlambdge0_i(i,:));
    
    L_i_m_n(i,:,:) = inv(squeeze(R_i_m_n(i,:,:)));
    LGE0_i_m_n(i,:,:) = L_i_m_n(i,find(squeeze(lambd_i_k(i,:)>=0)),:);
    LLE0_i_m_n(i,:,:) = L_i_m_n(i,find(squeeze(lambd_i_k(i,:)<=0)),:);
    for k=1:sn
        C_i_k_m_n(i,k,:,:) = squeeze(R_i_m_n(i,:,:))*reshape(D_k_m_n(k,:,:), [ sn sn ])*squeeze(L_i_m_n(i,:,:));
    end    
end
clear A_i_m_n R_i_m_n LAMBD_i_m_n L_i_m_n D_k_m_n ;
[sw1, sLGE0m, sLGE0n] = size(LGE0_i_m_n);
[sw1, sLLE0m, sLLE0n] = size(LLE0_i_m_n);

C_k1_k2_m_n = zeros(sn,sn,sn,sn);
for k1 = 1:sn
    for k2 = 1:sn
        C_k1_k2_m_n(k1,k2,:,:) = squeeze(C_i_k_m_n(1,k1,:,:))*squeeze(C_i_k_m_n(2,k2,:,:));
    end
end
C_k1k2n_m = reshape(permute(C_k1_k2_m_n, [ 1 2 4 3 ]), [ sn*sn*sn sn ]);

for g = 1:sg
    C_p1_p2_m_n{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn, sn);
end
for g = 1:sg
    for p1 = 2:(sp_g_i(g,1)-1)
        for p2 = 2:(sp_g_i(g,2)-1)
            C_p1_p2_m_n{g}(p1, p2, :, :) = eye(sn);
        end
    end
    
    for p1 = 1
        for p2 = 2:(sp_g_i(g,2)-1)
            for k1 = find(lambd_i_k(1,:)<=0)
                for k2 = 1:sn
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
            for k1 = find(lambd_i_k(1,:)>=0)
                for k2 = 1:sn
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
    
    for p1 = 2:(sp_g_i(g,1)-1)
        for p2 = 1
            for k1 = 1:sn
                for k2 = find(lambd_i_k(2,:)<=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
    for p1 = 2:(sp_g_i(g,1)-1)
        for p2 = sp_g_i(g,2)
            for k1 = 1:sn
                for k2 = find(lambd_i_k(2,:)>=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end

    for p1 = 1
        for p2 = 1
            for k1 = find(lambd_i_k(1,:)<=0)
                for k2 = find(lambd_i_k(2,:)<=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
    
    for p1 = 1
        for p2 = sp_g_i(g,2)
            for k1 = find(lambd_i_k(1,:)<=0)
                for k2 = find(lambd_i_k(2,:)>=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
    
    for p1 = sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
            for k1 = find(lambd_i_k(1,:)>=0)
                for k2 = find(lambd_i_k(2,:)>=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end

    for p1 = sp_g_i(g,1)
        for p2 = 1
            for k1 = find(lambd_i_k(1,:)>=0)
                for k2 = find(lambd_i_k(2,:)<=0)
                    C_p1_p2_m_n{g}(p1, p2, :, :) = C_p1_p2_m_n{g}(p1, p2, :, :) + C_k1_k2_m_n(k1,k2,:,:);
                end
            end
        end
    end
end
    

for g = 1:sg
    for i = 1:sdim
        hmax_g_i(g,i) = max(h_g_i_p{g,i});
        hmin_g_i(g,i) = min(h_g_i_p{g,i});
    end
end
hmax_i = max(hmax_g_i,[],1);
hmin_i = min(hmin_g_i,[],1);

deltat = min(hmin_i./transpose(max(abs(lambd_i_k),[],2)))*0.5;
Ntau = fix(T/deltat);
deltaNinput = fix((Ntau+1)/(Ninput-1));


%%  ”слови€ сопр€жени€ на границах
%   ћатрица свободных граничных условий
Afree_m_n = zeros(2,sn);    Afree_m_n(1,3) = 1; Afree_m_n(2,2) = 1;

%   ћатрица  условий слипани€
Aadhes_m_n = zeros(4,2*sn); 
Aadhes_m_n(1,1) = 1;    Aadhes_m_n(1,1+sn) = -1;
Aadhes_m_n(2,3) = 1;    Aadhes_m_n(2,3+sn) = -1;
Aadhes_m_n(3,4) = 1;    Aadhes_m_n(3,4+sn) = -1;
Aadhes_m_n(4,5) = 1;    Aadhes_m_n(4,5+sn) = -1;

Aadhesy_m_n = zeros(4,2*sn); 
Aadhesy_m_n(1,2) = 1;    Aadhesy_m_n(1,2+sn) = -1;
Aadhesy_m_n(2,3) = 1;    Aadhesy_m_n(2,3+sn) = -1;
Aadhesy_m_n(3,4) = 1;    Aadhesy_m_n(3,4+sn) = -1;
Aadhesy_m_n(4,5) = 1;    Aadhesy_m_n(4,5+sn) = -1;

A1adhes_m_n = zeros(7,2*sn); 
A1adhes_m_n(1,1) = 1;    A1adhes_m_n(1,1+sn) = -1;
A1adhes_m_n(2,3) = 1;    A1adhes_m_n(2,3+sn) = -1;
A1adhes_m_n(3,4) = 1;    A1adhes_m_n(3,4+sn) = -1;
A1adhes_m_n(4,5) = 1;    A1adhes_m_n(4,5+sn) = -1;
A1adhes_m_n(5,2) = 1;    
                         A1adhes_m_n(6,2+sn) = 1;
A1adhes_m_n(7,3) = 1;    
%                          A1adhes_m_n(8,3+sn) = 1;
rank(A1adhes_m_n)                         

Asliding_m_n = zeros(4,2*sn); 
Asliding_m_n(1,1)    = 1;   Asliding_m_n(1,1+sn) = -1;
Asliding_m_n(2,3)    = 1;    
Asliding_m_n(3,3+sn) = 1;    
Asliding_m_n(4,4) = 1;      Asliding_m_n(4,4+sn) = -1;

%%  ‘ормируем матрицы —Ћј” в точках границы    
R_g_p1_p2_m_n = cell(sg,max(sp_g_i(g,1)),max(sp_g_i(g,2)));

g = 1;
%   bound = 1   g = 1
    for p1 = 1
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 2   g = 1
    for p1 = 1:(sp_g_i(g,1))
        for p2 = sp_g_i(g,2)
            A = Afree_m_n;
            C = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C;
%   bound = 3   g = 1
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 1;   gpls = 2;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = A1adhes_m_n;
%             A = Asliding_m_n;
%             gmin = 1;   gpls = 2;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%             size(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 1
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 4;   gpls = 1;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
        
g = 2;
%   bound = 1   g = 2
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 1;   gpls = 2;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 2
    for p1 = 1:(sp_g_i(g,1))
        for p2 = sp_g_i(g,2)
            A = Afree_m_n;
            C = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
%   bound = 3   g = 2
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 2;   gpls = 3;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = Aadhes_m_n;
%             A = Asliding_m_n;
%             gmin = 2;   gpls = 3;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 2
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 5;   gpls = 2;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
        
g = 3;
%   bound = 1   g = 3
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 2;   gpls = 3;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 3
    for p1 = 1:sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
            A = Afree_m_n;
            C = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
%   bound = 3   g = 3
    for p1 = sp_g_i(g,1)
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 4
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 6;   gpls = 3;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
    
g = 4;
%   bound = 1   g = 4
    for p1 = 1
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 2   g = 4
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 4;   gpls = 1;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C;
%   bound = 3   g = 4
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 4;   gpls = 5;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = A1adhes_m_n;
%             A = Asliding_m_n;
%             gmin = 1;   gpls = 2;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%             size(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 4
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 7;   gpls = 4;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
        
g = 5;
%   bound = 1   g = 5
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 4;   gpls = 5;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 5
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 5;   gpls = 2;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C B;
%   bound = 3   g = 5
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 5;   gpls = 6;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = Aadhes_m_n;
%             A = Asliding_m_n;
%             gmin = 2;   gpls = 3;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 5
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 8;   gpls = 5;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;
        
g = 6;
%   bound = 1   g = 6
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 5;   gpls = 6;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 6
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 6;   gpls = 3;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C B;
%   bound = 3   g = 6
    for p1 = sp_g_i(g,1)
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 4
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
             A = Aadhesy_m_n;
%            A = Asliding_m_n;
            gmin = 9;   gpls = 6;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(2,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
    clear A C B;

    g = 7;
%   bound = 1   g = 7
    for p1 = 1
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 2   g = 7
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 7;   gpls = 4;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C;
%   bound = 3   g = 7
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 7;   gpls = 8;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = A1adhes_m_n;
%             A = Asliding_m_n;
%             gmin = 1;   gpls = 2;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%             size(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 7
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
        
g = 8;
%   bound = 1   g = 8
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 7;   gpls = 8;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 8
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 8;   gpls = 5;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C B;
%   bound = 3   g = 8
    for p1 = sp_g_i(g,1)
        for p2 = 2:(sp_g_i(g,2)-1)
             A = Aadhes_m_n;
%            A = Asliding_m_n;
            gmin = 8;   gpls = 9;
            C = zeros(sLGE0m+sLLE0m, 2*sn );
            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
        end
    end
%     for p1 = sp_g_i(g,1)
%         for p2 = sp_g_i(g,2)
% %             A = Aadhes_m_n;
%             A = Asliding_m_n;
%             gmin = 2;   gpls = 3;
%             C = zeros(sLGE0m+sLLE0m, 2*sn );
%             C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
% %             rank(cat(1,C,A))
%             R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%         end
%     end
    
    clear A C B;
%   bound = 4   g = 8
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
        
g = 9;
%   bound = 1   g = 9
%     for p1 = 1
%         for p2 = 1:sp_g_i(g,2)
%             A = Aadhes_m_n;
%             gmin = 8;   gpls = 9;
%             C = zeros(2*sn, 2*sn );
%             C(1:sn,1:sn) = squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :));
%             C((sn+1):(2*sn),(sn+1):(2*sn)) = squeeze(C_p1_p2_m_n{gpls}(1, p2, :, :));
%             B = cat(1,C,A);
%             R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
%         end
%     end
%     clear A C B;
%   bound = 2   g = 9
%    for p1 = 1:(sp_g_i(g,1))
%        for p2 = sp_g_i(g,2)
%             A = Aadhes_m_n;
%            A = Asliding_m_n;
%            gmin = 9;   gpls = 6;
%            C = zeros(sLGE0m+sLLE0m, 2*sn );
%            C(1:sLGE0m,1:sn) = squeeze(LGE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gmin}(p1, sp_g_i(gmin,2), :, :));
%            C((sLGE0m+1):(sLGE0m+sLLE0m),(sn+1):(2*sn)) = squeeze(LLE0_i_m_n(1,:,:))*squeeze(C_p1_p2_m_n{gpls}(p1, 1, :, :));
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(squeeze(C_p1_p2_m_n{gmin}(sp_g_i(gmin,1), p2, :, :)))
%             rank(cat(1,C,A))
%            R_g_p1_p2_m_n{g,p1,p2} = inv(cat(1,C,A));
%        end
%    end
%    clear A C B;
%   bound = 3   g = 9
    for p1 = sp_g_i(g,1)
        for p2 = 1:(sp_g_i(g,2)-1)
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;
%   bound = 4
    for p1 = 2:(sp_g_i(g,1))
        for p2 = 1
            %   ѕрозрачные граничные услови€
            A = eye(sn) - squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            C = squeeze(C_p1_p2_m_n{g}(p1,p2, :, :));
            B = cat(1,C,A);
            R_g_p1_p2_m_n{g,p1,p2} = inv(transpose(B)*B)*transpose(B);
        end
    end
    clear A C B;

%%  ‘ормируем базисные полиномы
H_g_i_pw_k_p = cell(sg,sdim);
for g = 1:sg
    for i = 1:sdim
        H_g_i_pw_k_p{g,i} = zeros(sp_g_i(g,i),sn,sp_g_i(g,i));
    end
end

for g = 1:sg
    for i = 1:sdim
%         H_pw_k_p = zeros(sp_g_i(g,i),sn,sp_g_i(g,i));
%         h_p = h_g_i_p{g,i};
        for pw = 1
            for k = find(lambd_i_k(i,:)<=0)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(1,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_i_k(i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw+1) = polyval(pbase_l_k(2,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_i_k(i,k)*deltat));
            end
        end

        for pw = 2:(sp_g_i(g,i)-1)
            for k = find(lambd_i_k(i,:)<=0)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(1,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_i_k(i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw+1) = polyval(pbase_l_k(2,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_i_k(i,k)*deltat));
            end
            for k = find(lambd_i_k(i,:)>=0)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(2,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_i_k(i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw-1) = polyval(pbase_l_k(1,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_i_k(i,k)*deltat));
            end
        end

        for pw = sp_g_i(g,i)
            for k = find(lambd_i_k(i,:)>=0)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(2,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_i_k(i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw-1) = polyval(pbase_l_k(1,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_i_k(i,k)*deltat));
            end
        end
    end
end

for g = 1:sg
    for i = 1:sdim
        H_g_i_pwk_p{g,i} = sparse(reshape(H_g_i_pw_k_p{g,i}, [ sp_g_i(g,i)*sn sp_g_i(g,i)]));
    end
end

% for g = 1:sg
%     for i = 1:sdim
%         H_g_i_pwp_k{g,i} = sparse(reshape(permute(H_g_i_pw_k_p{g,i}, [ 1 3 2 ]), [ sp_g_i(g,i)*sp_g_i(g,i) sn ]));
%     end
% end
clear H_g_i_pw_k_p;

C_i_k_mn = reshape(C_i_k_m_n, [ sdim sn sn*sn ]);
clear C_i_k_m_n;

% for i = 1:sdim
%     indk_i_k(i, 1:sn) = 1:sn;
% end
% 
% indspm_g_i_m = zeros(sg,sdim,sn);
% for g = 1:sg
%     for i = 1:sdim
%         for m = 1:sn
%             indspm_g_i_m(g,i,m) = sub2ind([sp_g_i(g,i),sn],sp_g_i(g,i),m);
%             ind1m_g_i_m(g,i,m) = sub2ind([sp_g_i(g,i),sn],1,m);
%             
%         end
%     end
% end
% 
% for g = 1:sg
%     for i = 1:sdim
%         HC_g_i_pwm_pn{g,i} = sparse(reshape(...
%                                         permute(...
%                                             reshape(H_g_i_pwp_k{g,i}*squeeze(C_i_k_mn(i,:,:)), ...
%                                             [ sp_g_i(g,i) sp_g_i(g,i) sn sn ]), ...
%                                         [ 1 3 2 4 ]), ...
%                                     [ sp_g_i(g,i)*sn sp_g_i(g,i)*sn ]));
%                                 
% %         HCLGE0_g_i_pwm_pn{g,i} = sparse(reshape(...
% %                                         permute(...
% %                                             reshape(H_g_i_pwp_k{g,i}(:,indk_i_k(i,lambd_i_k(i,:)>=0))*squeeze(C_i_k_mn(i,:,:)), ...
% %                                             [ sp_g_i(g,i) sp_g_i(g,i) sn sn ]), ...
% %                                         [ 1 3 2 4 ]), ...
% %                                     [ sp_g_i(g,i)*sn sp_g_i(g,i)*sn ]));
% %                                 
% %         HCLGE0_g_i_spm_pn{g,i} = sparse(HCLGE0_g_i_pwm_pn{g,i}(squeeze(indspm_g_i_m(g,i,:)),:));
% %                                 
% %         HCLLE0_g_i_pwm_pn{g,i} = sparse(reshape(...
% %                                         permute(...
% %                                             reshape(H_g_i_pwp_k{g,i}(:,indk_i_k(i,lambd_i_k(i,:)<=0))*squeeze(C_i_k_mn(i,:,:)), ...
% %                                             [ sp_g_i(g,i) sp_g_i(g,i) sn sn ]), ...
% %                                         [ 1 3 2 4 ]), ...
% %                                     [ sp_g_i(g,i)*sn sp_g_i(g,i)*sn ]));
% %                                 
% %         HCLLE0_g_i_1m_pn{g,i} = sparse(HCLLE0_g_i_pwm_pn{g,i}(squeeze(ind1m_g_i_m(g,i,:)),:));
%     end
% end
% clear H_g_i_pwp_k;

%%  Ќачальные данные u0_p1_p2_n
for g = 1:sg
    u0_g_p1_p2_n{g} = zeros( sp_g_i(g,1),sp_g_i(g,2),sn );
end

g = 2;
    u0_p1_p2_n = zeros( sp_g_i(g,1),sp_g_i(g,2),sn );
    amplitudeu0 = 0;
    xstar_i = [0 0];
    diam = 0.1;
    centerline = 0.3;
    x_p1_p2_i = zeros(sp_g_i(g,1),sp_g_i(g,2),sdim);
    for p2 = 1:sp_g_i(g,2)
        x_p1_p2_i(:,p2,1) = reshape(x_p{g,1}, [ sp_g_i(g,1) 1 1 ]);
    end
    for p1 = 1:sp_g_i(g,1)
        x_p1_p2_i(p1,:,2) = reshape(x_p{g,2}, [ 1 sp_g_i(g,2) 1 ]);
    end

    radius_p1_p2_i = zeros(sp_g_i(g,1),sp_g_i(g,2),sdim);
    normradius_p1_p2 = zeros(sp_g_i(g,1),sp_g_i(g,2));
for p1 = 1:sp_g_i(g,1)
    for p2 = 1:sp_g_i(g,2)
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
    clear indp1 indp2;
    
    u0_g_p1_p2_n{g} = u0_p1_p2_n;
% for g = 1:sg
%     u0_g_p_n{g} = zeros( sp_g_i(g,:),sn );
% end
% 
% g = 2;
% amplitudeu0 = 10;
% xstar_p_i = zeros(sp_g_i(g,:),sdim);
% diam = 0.1;
% centerline = 0.3;
% 
% 
% for i = 1:sdim
%     x_p_i(1:sp_g_i(g,i),i) = x_p{g,i};
% end
% radius_p_i = x_p_i - xstar_p_i;
% normradius_p = (sum(radius_p_i.^2, sdim+1).^(0.5));
% 
% ind_p = find(abs(normradius_p - centerline) <= diam);
% lengthindp = length(ind_p);
% u0_p_n = zeros(sp_g_i(g,:),sn);
% for point = 1:lengthindp
%     u0_p_n(ind2sub(sp_g_i(g,:),ind_p(point),4:5)) = amplitudeu0*...
%        (radius_p_i(ind2sub(sp_g_i(g,:),ind_p(point),4:5),:)/normradius_p(ind2sub(sp_g_i(g,:),ind_p(point))))...
%        *(1 + cos((pi / diam)*(normradius_p(ind2sub(sp_g_i(g,:),ind_p(point)))-centerline)));
% end
% clear ind_p;
% 
%%   ƒелаем шаги по времени
toc;
tic;

sjt = Ninput;
T_jt = 0:deltat*deltaNinput:(sjt-1)*deltat*deltaNinput;
indk1 = 1:sn;  indk2 = 1:sn;

for g = sg:-1:1
    amplitude_p1w_sp2{g} = zeros(sp_g_i(g,1),1);
    g = 2;
    amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2-1, 1) = 100;    
        amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2, 1) = 100;
            amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2+1, 1) = 100;

    uprev_g_p1_p2_n = u0_g_p1_p2_n;
    
    u_g_p1_p2_n_jt{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn, sjt);
    u_g_p1_p2_n_jt{g}(:,:,:,1) = u0_g_p1_p2_n{g};

%     unext_p1_p2_n = zeros(sp_g_i(g,1), sp_g_i(g,2),sn);
%     W_m_p1w = zeros(sn,sp_g_i(g,1));
end

for g = sg:-1:1
    unext_g_p1_p2_n{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn);
end
for jt=2:sjt
    for point = 2:deltaNinput+1
        for g = 1:sg
            %   ¬ычисл€ем значени€ во всех узлах без учета условий на
            %   внутренних\внешних границах
            uprev_g_p1_p2n = reshape( uprev_g_p1_p2_n{g}, [ sp_g_i(g,1) sp_g_i(g,2)*sn ]);
            W_p1w_k1_p2_n= reshape(H_g_i_pwk_p{g,1}*uprev_g_p1_p2n, [ sp_g_i(g,1) sn sp_g_i(g,2) sn ]);
            W_p2w_k2_p1w_k1_n = reshape(H_g_i_pwk_p{g,2}*reshape(permute(W_p1w_k1_p2_n, [ 3 1 2 4 ]), ...
                                                      [ sp_g_i(g,2) sp_g_i(g,1)*sn*sn ]), ...
                                                            [ sp_g_i(g,2) sn sp_g_i(g,1) sn sn ]);
            unext_g_p1_p2_n{g} = reshape(...
                reshape(permute(W_p2w_k2_p1w_k1_n, [ 3 1 4 2 5 ]), [ sp_g_i(g,1)*sp_g_i(g,2) sn*sn*sn ])...
                *C_k1k2n_m, [ sp_g_i(g,1) sp_g_i(g,2) sn ]);                                            
        end

        g = 1;
        %   bound = 2   g = 1
        for p1 = 1:(sp_g_i(g,1))
            for p2 = sp_g_i(g,2)
                r1_m =squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{g}(p1,p2,:));
                t = (jt-1)*deltat*deltaNinput + point*deltat;
                r2_m(1,1) = 0; ...
                    r2_m(2,1) = -(g==2)*transpose(amplitude_p1w_sp2{g}(p1,1)*sin(omega*t));
                rx_m = cat(1,r1_m,r2_m); 
                unext_g_p1_p2_n{g}(p1,p2,:) = R_g_p1_p2_m_n{g,p1,p2}*rx_m;
            end
        end
        %   bound = 3   g = 1
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 1;   gpls = 2;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end
        %   bound = 4   g = 1
        for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 4;   gpls = 1;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 1;   gpls = 2;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m(1:4,1) = 0; 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 2;
        %   bound = 2   g = 2
        for p1 = 1:(sp_g_i(g,1))
            for p2 = sp_g_i(g,2)
                r1_m =squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{g}(p1,p2,:));
                t = (jt-1)*deltat*deltaNinput + point*deltat;
                r2_m(1,1) = 0; ...
                    r2_m(2,1) = -(g==2)*transpose(amplitude_p1w_sp2{g}(p1,1)*sin(omega*t));
                rx_m = cat(1,r1_m,r2_m); 
                unext_g_p1_p2_n{g}(p1,p2,:) = R_g_p1_p2_m_n{g,p1,p2}*rx_m;
            end
        end
    %   bound = 3     g = 2
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 2;   gpls = 3;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m = zeros(4,1); 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end
        
        for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 5;   gpls = 2;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 2;   gpls = 3;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m = zeros(4,1); 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 3;
        %   bound = 2     g = 3
        for p1 = 1:sp_g_i(g,1)
            for p2 = sp_g_i(g,2)
                r1_m =squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{g}(p1,p2,:));
                t = (jt-1)*deltat*deltaNinput + point*deltat;
                r2_m(1,1) = 0; ...
                    r2_m(2,1) = -(g==2)*transpose(amplitude_p1w_sp2{g}(p1,1)*sin(omega*t));
                rx_m = cat(1,r1_m,r2_m); 
                unext_g_p1_p2_n{g}(p1,p2,:) = R_g_p1_p2_m_n{g,p1,p2}*rx_m;
            end
        end
  
         for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 6;   gpls = 3;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
        
 
        g = 4;
        %   bound = 3   g = 1
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 4;   gpls = 5;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end
        %   bound = 4   g = 4
        for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 7;   gpls = 4;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 1;   gpls = 2;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m(1:4,1) = 0; 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 5;
    %   bound = 3     g = 5
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 5;   gpls = 6;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m = zeros(4,1); 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end
        
         for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 8;   gpls = 5;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 2;   gpls = 3;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m = zeros(4,1); 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 6;
        for p1 = 2:(sp_g_i(g,1)-1)
            for p2 = 1
                gmin = 9;   gpls = 6;
                r11_m = squeeze(LGE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(2,:,:))*squeeze(unext_g_p1_p2_n{gpls}(p1,1,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(p1,1,:) = W(sn+1:sn+sn,1);
            end
        end
    
 
        g = 7;
        %   bound = 3   g = 7
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 1;   gpls = 2;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m(1:4,1) = 0; 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end

%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 1;   gpls = 2;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m(1:4,1) = 0; 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 8;
    %   bound = 3     g = 8
        for p1 = sp_g_i(g,1)
            for p2 = 2:(sp_g_i(g,2)-1)
                gmin = 8;   gpls = 9;
                r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
                r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
                r21_m = zeros(4,1); 
                rxx_m = cat(1,r11_m,r12_m,r21_m);
                W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
                unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
            end
        end
%         for p1 = sp_g_i(g,1)
%             for p2 = sp_g_i(g,2)
%                 gmin = 2;   gpls = 3;
%                 r11_m = squeeze(LGE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gmin}(p1,p2,:));
%                 r12_m = squeeze(LLE0_i_m_n(1,:,:))*squeeze(unext_g_p1_p2_n{gpls}(1,p2,:));
%                 r21_m = zeros(4,1); 
%                 rxx_m = cat(1,r11_m,r12_m,r21_m);
%                 W = R_g_p1_p2_m_n{g,p1,p2}*rxx_m;
%                 unext_g_p1_p2_n{gmin}(p1,p2,:) = W(1:sn,1); unext_g_p1_p2_n{gpls}(1,p2,:) = W(sn+1:sn+sn,1);
%             end
%         end
        
        g = 9;
        
        uprev_g_p1_p2_n = unext_g_p1_p2_n;
    end
%     for g = 1:sg
%         u_g_p1_p2_n_jt{g}(:,:,:,jt) = unext_g_p1_p2_n{g};
%     end
u_g_p1_p2_n_jt{1}(:,:,:,jt) = unext_g_p1_p2_n{1};
u_g_p1_p2_n_jt{2}(:,:,:,jt) = unext_g_p1_p2_n{2};
u_g_p1_p2_n_jt{3}(:,:,:,jt) = unext_g_p1_p2_n{3};
u_g_p1_p2_n_jt{4}(:,:,:,jt) = unext_g_p1_p2_n{4};
u_g_p1_p2_n_jt{5}(:,:,:,jt) = unext_g_p1_p2_n{5};
u_g_p1_p2_n_jt{6}(:,:,:,jt) = unext_g_p1_p2_n{6};
u_g_p1_p2_n_jt{7}(:,:,:,jt) = unext_g_p1_p2_n{7};
u_g_p1_p2_n_jt{8}(:,:,:,jt) = unext_g_p1_p2_n{8};
u_g_p1_p2_n_jt{9}(:,:,:,jt) = unext_g_p1_p2_n{9};
end


for g = 1:sg
    K_g_p1_p2_jt{g} = zeros(sp_g_i(gmin,1),sp_g_i(gmin,2),sjt);
end

for g = 1:sg
    K_g_p1_p2_jt{g} = (squeeze(u_g_p1_p2_n_jt{g}(:,:,4,:)).^2 + squeeze(u_g_p1_p2_n_jt{g}(:,:,5,:)).^2).^(0.5);
end

toc;

%%  √рафический вывод результатов
x_p1 = cat(2,x_p{1,1},x_p{2,1},x_p{3,1});   
x_p2 = x_p{1,2};
K_p1_p2_jt = cat(1,K_g_p1_p2_jt{1},K_g_p1_p2_jt{2},K_g_p1_p2_jt{3})*10000;
[sp1,sp2,sj] = size(K_p1_p2_jt);

x1min = min(x_p1);                       x1max = max(x_p1);
x2min = min(x_p2);                       x2max = max(x_p2);

Kmin = min(min(min(K_p1_p2_jt, [], 1), [], 2), [], 3);   Kmax = max(max(max(K_p1_p2_jt, [], 1), [], 2), [], 3);

[X2, X1] = meshgrid(x_p2, x_p1);

[per, sit] = size(T_jt);    clear per;

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
            R(1, :, :) = X2;
            R(2, :, :) = X1;
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

% %%  √рафический вывод результатов
% for g = 1:sg
%     for i = 1:sdim
%         xmin_g_i(g,i) = min(x_p{g,i});  xmax_g_i(g,i) = max(x_p{g,i});
%     end
% end
% xmin_i = min(xmin_g_i,[],1);    xmax_i = max(xmax_g_i,[],1);
% 
% for g = 1:sg
%     Kmin_g(g,1) = min(min(min(K_g_p1_p2_jt{g}, [], 1), [], 2), [], 3);  Kmax_g(g,1) = max(max(max(K_g_p1_p2_jt{g}, [], 1), [], 2), [], 3);
% end
% Kmin = min(Kmin_g); Kmax = max(Kmax_g);
% 
% [per, sit] = size(T_jt);    clear per;
% 
% for it=1:1:sit
%     for g = 1:sg
%         surf(x_p{g,2},x_p{g,1},reshape(K_g_p1_p2_jt{g}(:,:,it), [ sp_g_i(g,1) sp_g_i(g,2) ]));
%         axis([ xmin_i(2) xmax_i(2) xmin_i(1) xmax_i(1) Kmin Kmax]);
%         hold on;
%     end
%     hold off;
%     pause(0.01);
%     F(it) = getframe; %#ok<AGROW>
% end
%         
% movie(F, 3, 5);
        
% function wbound_g_p1_p2_n = Bound_external(unext_g_i_bound_p1_p2_n, i, bound, cond)
% switch i
%     case 1
%         switch bound
%             case 1
%                 switch cond
%                     case 1
%                     case 2
%                 end
%             case 2
%                 switch cond
%                     case 1
%                     case 2
%                 end
%         end
%     case 2
%         switch bound
%             case 1
%                 switch cond
%                     case 1
%                     case 2
%                 end
%             case 2
%                 switch cond
%                     case 1
%                     case 2
%                         W_m_p1w = transpose(squeeze(unext_g_i_bound_p1_p2_n{g,i,bound}));
%                         r1_m_p1w = squeeze(LGE0_i_m_n(2,:,:))*W_m_p1w;
%                         t = (jt-1)*deltat*deltaNinput + point*deltat;
%                         r2_m_p1w(1,1:sp_g_i(g,1)) = 0; ...
%                             r2_m_p1w(2,1:sp_g_i(g,1)) = -(g==2)*transpose(amplitude_p1w_sp2(1:sp_g_i(g,1),1)*sin(omega*t));
%                         r_m_p1w = cat(1,r1_m_p1w,r2_m_p1w); 
%                         wnext_g_i_bound_p1_p2_n = reshape((invBext_i_cond_m_n{2,1}*r_m_p1w), [ sp_g_i(g,1) 1 sn ]);
%                 end
%         end
% end
%         
% end
% 
% 
% function [wbound_gmin_p1_p2_n, wbound_gpls_p1_p2_n] = Bound_internal(unext_g_i_bound_p1_p2_n, i, gmin, gpls, cond)
% switch i
%     case 1
%         switch cond
%             case 1
%                 W_m_p2w = transpose(squeeze(unext_g_i_bound_p1_p2_n{gmin,i,2}));
%                 r11_m_p2w = squeeze(LGE0_i_m_n(1,:,:))*W_m_p1w;
%                 W_m_p2w = transpose(squeeze(unext_g_i_bound_p1_p2_n{gpls,i,1}));
%                 r12_m_p2w = squeeze(LLE0_i_m_n(1,:,:))*W_m_p1w;
%                 
%                 r2_m_p2w(1:4,1:sp_g_i(g,2)) = 0; 
%                 r_m_p2w = cat(1,r11_m_p2w,r12_m_p2w,r2_m_p2w);
%                 
%                 wnext_p1_p2_n = reshape((invBint_i_cond_m_n{i,cond}*r_m_p2w), [ 1 2*sp_g_i(g,2) sn ]);
%                 wbound_gmin_p1_p2_n = wnext_p1_p2_n(:,1:sp_g_i(g,2),:);
%                 wbound_gpls_p1_p2_n = wnext_p1_p2_n(:,sp_g_i(g,2)+1:2*sp_g_i(g,2),:);
%             case 2
%                 W_m_p2w = transpose(squeeze(unext_g_i_bound_p1_p2_n{gmin,i,2}));
%                 r11_m_p2w = squeeze(LGE0_i_m_n(1,:,:))*W_m_p1w;
%                 W_m_p2w = transpose(squeeze(unext_g_i_bound_p1_p2_n{gpls,i,1}));
%                 r12_m_p2w = squeeze(LLE0_i_m_n(1,:,:))*W_m_p1w;
%                 
%                 r2_m_p2w(1:4,1:sp_g_i(g,2)) = 0; 
%                 r_m_p2w = cat(1,r11_m_p2w,r12_m_p2w,r2_m_p2w);
%                 
%                 wnext_p1_p2_n = reshape((invBint_i_cond_m_n{i,cond}*r_m_p2w), [ 1 2*sp_g_i(g,2) sn ]);
%                 wbound_gmin_p1_p2_n = wnext_p1_p2_n(:,1:sp_g_i(g,2),:);
%                 wbound_gpls_p1_p2_n = wnext_p1_p2_n(:,sp_g_i(g,2)+1:2*sp_g_i(g,2),:);
%         end
%     case 2
%         switch cond
%             case 1
%             case 2
%         end
% end
%         
% end

%end


