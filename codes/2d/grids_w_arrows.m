% Этот файл получен из global_init2_1_3
% Адаптирован для решения системы уравнений динамической упругости
% многомерной модификацией сеточно-характеристического метода, 
% использующей фундаментальное решение оператора задачи.

function grid_charact_2_1_1_1_4_m()
format long;
%% Блок для описания начальных данных:
%Засекаем время работы программы:
tic;
profile on
sdim = 2;                                                                   %Размерность по пространственным переменым
sn = 5;
sg = 3;

%   Задаем параметры задачи
for g = 1:sg
    la_g(g) = 2; mu_g(g) = 1; ro_g(g) = 1;                        % параметры среды
end

A_g_i_m_n = zeros(sg,sdim,sn,sn);
for g = 1:sg
A_g_i_m_n(g,1,:,:) = ...
[ 0         0   0       -(la_g(g)+2*mu_g(g))  0   ;...
  0         0   0       -la_g(g)         0   ;...
  0         0   0       0           -mu_g(g) ;...
  -1/ro_g(g)     0   0       0           0   ;...
  0         0   -1/ro_g(g)   0           0   ];

A_g_i_m_n(g,2,:,:) = ...
[ 0 0       0       0   -la_g(g)         ;...
  0 0       0       0   -(la_g(g)+2*mu_g(g))  ;...
  0 0       0       -mu_g(g) 0           ;...
  0 0       -1/ro_g(g)   0   0           ;...
  0 -1/ro_g(g)   0       0   0           ];
end
                                                                      % количество регионов

omega =1;                                                                  %частота действующей силы
% %Временные параметры:
T = 5.1;
Ninput= 61;

%   Задаем область, в которой ищем решение
dmin_i = zeros(1,sdim); dmax_i = zeros(1,sdim);
dmin_i(1) = -1;  dmax_i(1) = 1;
dmin_i(2) = -1;  dmax_i(2) = 1;

dmin_g_i = zeros(sg,sdim); dmax_g_i = zeros(sg,sdim);
dmin_g_i(1,1) = dmin_i(1); dmax_g_i(1,1) = dmin_i(1)/2;
    dmin_g_i(2,1) = dmin_i(1)/2; dmax_g_i(2,1) = dmax_i(1)/2;
        dmin_g_i(3,1) = dmax_i(1)/2; dmax_g_i(3,1) = dmax_i(1);             % Левая и правая границы области
dmin_g_i(1,2) = dmin_i(2); dmax_g_i(1,2) = dmax_i(2);
    dmin_g_i(2,2) = dmin_i(2); dmax_g_i(2,2) = dmax_i(2);
        dmin_g_i(3,2) = dmin_i(2); dmax_g_i(3,2) = dmax_i(2);
    
%   Задаем порядок интерполяционных полиномов по пространству
Kmax = 1;  

%   Задаем сетку по пространственным переменным x_ii1 и x_ii2
sp_i = zeros(1,sdim);   sp_i(1) = 2*12+1;  sp_i(2) = 2*12+1;              % Количество узлов сетки-всегда нечетное
sp_g_i = zeros(sg,sdim);
sp_g_i(1,1) = (sp_i(1)-1)/4 + 1;
    sp_g_i(2,1) = 2*(sp_i(1)-1)/4 + 1;
        sp_g_i(3,1) = (sp_i(1)-1)/4 + 1;
sp_g_i(1,2) = sp_i(2);
    sp_g_i(2,2) = sp_i(2);
        sp_g_i(3,2) = sp_i(2);

        
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

for g = 1:sg
    for i = 1:sdim
        hmax_g_i(g,i) = max(h_g_i_p{g,i});
        hmin_g_i(g,i) = min(h_g_i_p{g,i});
    end
end

% Массивы C_g_i_k_m_n
D_k_m_n = zeros(sn,sn,sn);
for k=1:sn
    D_k_m_n(k,k,k) = 1;
end
C_g_i_k_m_n = zeros(sg,sdim,sn,sn,sn);

for g = 1:sg
    for i = 1:sdim
    [R_g_i_m_n(g,i,:,:),LAMBD_g_i_m_n(g,i,:,:)] = eig(squeeze(A_g_i_m_n(g,i,:,:)));
    lambd_g_i_k(g,i,:) = diag(squeeze(LAMBD_g_i_m_n(g,i,:,:)));
    end
end

tollambd = 10^(-10);
for g = 1:sg
    for i = 1:sdim
        L_g_i_m_n(g,i,:,:) = inv(squeeze(R_g_i_m_n(g,i,:,:)));
        LGT0_g_i_m_n(g,i,:,:) = L_g_i_m_n(g,i,find(squeeze(lambd_g_i_k(g,i,:))>  tollambd),:);
        LGE0_g_i_m_n(g,i,:,:) = L_g_i_m_n(g,i,find(squeeze(lambd_g_i_k(g,i,:))>=-tollambd),:);
        LLE0_g_i_m_n(g,i,:,:) = L_g_i_m_n(g,i,find(squeeze(lambd_g_i_k(g,i,:))<= tollambd),:);
        LLT0_g_i_m_n(g,i,:,:) = L_g_i_m_n(g,i,find(squeeze(lambd_g_i_k(g,i,:))< -tollambd),:);
        for k=1:sn
            C_g_i_k_m_n(g,i,k,:,:) = squeeze(R_g_i_m_n(g,i,:,:))...
                                    *reshape(D_k_m_n(k,:,:), [ sn sn ])...
                                    *squeeze(L_g_i_m_n(g,i,:,:));
        end
    end
end

hmax_i = max(hmax_g_i,[],1);
hmin_i = min(hmin_g_i,[],1);
lwave = max(max(max(abs(lambd_g_i_k))))*2*pi/omega;

delta1 = min( hmin_i./max(max(max(abs(lambd_g_i_k)))));
delta2 = (2*pi/omega)/10;
deltat = min(delta2,delta1);
Ntau = fix(T/deltat)
deltaNinput = fix((Ntau+1)/(Ninput-1))

sc = prod(1:sdim);
ind_i_c = zeros(sdim,sc);    ind_i_c(:,1) = [1;2];  ind_i_c(:,2) = [2;1];

C_g_c_k1_k2_m_n = cell(sg,sc);
for g = 1:sg
    for c = 1:sc
        C_g_c_k1_k2_m_n(g,c) = {zeros(sn,sn,sn,sn)};
    end
end

for g = 1:sg
    for c = 1:sc
        for k1 = 1:sn
            for k2 = 1:sn
                k_i = [k1;k2];
                C_g_c_k1_k2_m_n{g,c}(k1,k2,:,:) = squeeze(C_g_i_k_m_n(g,ind_i_c(1,c),k_i(ind_i_c(1,c)),:,:))...
                                                 *squeeze(C_g_i_k_m_n(g,ind_i_c(2,c),k_i(ind_i_c(2,c)),:,:));
            end
        end
    end
end
C_g_c_k1k2n_m = cell(sg,sc);
for g = 1:sg
    for c = 1:sc
        C_g_c_k1k2n_m(g,c) = {sparse(reshape(permute(C_g_c_k1_k2_m_n{g,c}, [ 1 2 4 3 ]), [ sn*sn*sn sn ]))};
    end
end

%%  Формируем вспомагательные массивы ячеек    

C_g_c_p1_p2_m_n = internal_characteristic;        %Матрицы соответствующие внутреним характеристикам

B_g_c_bound_ind = boundary_conditions;            %Условия сопряжения на границах (внутренних/внешних)
    
H_g_i_pw_k_p = Basic_polynoms();                  %  Формируем базисные полиномы  

[R_g_c_p1_p2,Rd_g_c_p1_p2] = SLAE_boundary();          %  Формируем матрицы СЛАУ в точках границы


%%  Начальные данные u0_p1_p2_n
for g = 1:sg
    u0_g_p1_p2_n{g} = zeros( sp_g_i(g,1),sp_g_i(g,2),sn );
end

g = 2;
    u0_p1_p2_n = zeros( sp_g_i(g,1),sp_g_i(g,2),sn );
    amplitudeu0 = 100;
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

toc;
tic;
%%   Делаем шаги по времени

sjt = Ninput;
T_jt = 0:deltat*deltaNinput:(sjt-1)*deltat*deltaNinput;
indk1 = 1:sn;  indk2 = 1:sn;

for g = sg:-1:1
    amplitude_p1w_sp2{g} = zeros(sp_g_i(g,1),1);
    g = 2;
    amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2-1, 1) = 0;    
        amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2, 1) = 0;
            amplitude_p1w_sp2{g}((sp_g_i(g,1)+1)/2+1, 1) = 0;

    uprev_g_p1_p2_n = u0_g_p1_p2_n;
    
    u_g_p1_p2_n_jt{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn, sjt);
    u_g_p1_p2_n_jt{g}(:,:,:,1) = u0_g_p1_p2_n{g};

%     unext_p1_p2_n = zeros(sp_g_i(g,1), sp_g_i(g,2),sn);
%     W_m_p1w = zeros(sn,sp_g_i(g,1));
end

unext_g_c_p1_p2_n = cell(sg,sc);
unext_g_p1_p2_n = cell(sg);
for g = 1:sg
    unext_g_p1_p2_n{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn);
    for c = 1:sc
        unext_g_c_p1_p2_n{g,c} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn);
    end
end
for jt=2:sjt
    for point = 2:deltaNinput+1
        for g = 1:sg
            for c = 1:sc
            %   Вычисляем значения во всех узлах без учета условий на
            %   внутренних\внешних границах
            uprev_g_p1_p2n = reshape( uprev_g_p1_p2_n{g}, [ sp_g_i(g,1) sp_g_i(g,2)*sn ]);
            W_p1w_k1_p2_n= reshape(H_g_i_pwk_p{g,1}*uprev_g_p1_p2n, [ sp_g_i(g,1) sn sp_g_i(g,2) sn ]);
            W_p2w_k2_p1w_k1_n = reshape(H_g_i_pwk_p{g,2}*reshape(permute(W_p1w_k1_p2_n, [ 3 1 2 4 ]), ...
                                                      [ sp_g_i(g,2) sp_g_i(g,1)*sn*sn ]), ...
                                                            [ sp_g_i(g,2) sn sp_g_i(g,1) sn sn ]);
            unext_g_c_p1_p2_n{g,c} = reshape(...
                reshape(permute(W_p2w_k2_p1w_k1_n, [ 3 1 4 2 5 ]), [ sp_g_i(g,1)*sp_g_i(g,2) sn*sn*sn ])...
                *C_g_c_k1k2n_m{g,c}, [ sp_g_i(g,1) sp_g_i(g,2) sn ]);
            end
        end

for c = 1:sc        
g = 1;
%   bound = 1   bound = 4   g = 1
    p1 = 1; p2 = 1;
        gint = g;   gext = B_g_c_bound_ind{g,c,1,2};
        rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))];
        W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
        unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W;
%   bound = 1   g = 1
    p1 = 1; p2 = 2:(sp_g_i(g,2)-1);
            %   Прозрачные граничные условия 
        gint = g;   gext = B_g_c_bound_ind{g,c,1,2};
        rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
        W = R_g_c_p1_p2{g,c}{p1,2}*rxx_m ;
        unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[1 sp_g_i(g,2)-2 sn ]);
%   bound = 1   bound = 2   g = 1
    p1 = 1; p2 = sp_g_i(g,2);
            %   Прозрачные граничные условия + Условия свободной границы
        gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
        rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))];
        W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
        unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W;
%   bound = 2   g = 1
    p1 = 2:(sp_g_i(g,1)-1); p2 = sp_g_i(g,2);
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);
%   bound = 2   bound = 3   g = 1
    p1 = sp_g_i(g,1);   p2 = sp_g_i(g,2);
            %   Условия свободной границы + Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:));...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W(1:sn,1); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = W(sn+1:sn+sn,1);
%   bound = 3   g = 1
    p1 = sp_g_i(g,1);   p2 = 2:(sp_g_i(g,2)-1);
            %   Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))';...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))'];
            W = R_g_c_p1_p2{g,c}{p1,2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W(1:sn,:)',[1 sp_g_i(g,2)-2 sn ]); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = reshape(W(sn+1:sn+sn,:)',[1 sp_g_i(g,2)-2 sn ]);
%   bound = 3   bound = 4   g = 1
    p1 = sp_g_i(g,1);   p2 = 1;
            %   Условия проскальзывания/слипания + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:));...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W(1:sn,1); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = W(sn+1:sn+sn,1);
%   bound = 4   g = 1
    p1 = 2:(sp_g_i(g,1)-1); p2 = 1;
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);
       
g = 2;
%   bound = 2   g = 2
    p1 = 2:(sp_g_i(g,1)-1); p2 = sp_g_i(g,2);
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
            t = (jt-1)*deltat*deltaNinput + point*deltat;
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m + Rd_g_c_p1_p2{g,c}{2,p2}*d_g_p1_p2{g}{2,p2}*amplitude_p1w_sp2{g}(p1,1)'*sin(omega*t);
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);
%   bound = 2   bound = 3   g = 2
    p1 = sp_g_i(g,1);   p2 = sp_g_i(g,2);
            %   Условия свободной границы + Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:));...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W(1:sn,1); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = W(sn+1:sn+sn,1);
%   bound = 3   g = 2
    p1 = sp_g_i(g,1);   p2 = 2:(sp_g_i(g,2)-1);
            %   Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))';...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))'];
            W = R_g_c_p1_p2{g,c}{p1,2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W(1:sn,:)',[1 sp_g_i(g,2)-2 sn ]); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = reshape(W(sn+1:sn+sn,:)',[1 sp_g_i(g,2)-2 sn ]);
%   bound = 3   bound = 4   g = 2
    p1 = sp_g_i(g,1);   p2 = 1;
            %   Условия проскальзывания/слипания + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:));...
                     squeeze(unext_g_c_p1_p2_n{gext,c}(1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W(1:sn,1); unext_g_c_p1_p2_n{gext,c}(1,p2,:) = W(sn+1:sn+sn,1);
%   bound = 4   g = 2
    p1 = 2:(sp_g_i(g,1)-1); p2 = 1;
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m ;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);
        
g = 3;
%   bound = 2   g = 3'
    p1 = 2:(sp_g_i(g,1)-1); p2 = sp_g_i(g,2);
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))]';
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);
%   bound = 2   bound = 3   g = 3
    p1 = sp_g_i(g,1);   p2 = sp_g_i(g,2);
            %   Условия свободной границы + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W;
%   bound = 3   g = 3
    p1 = sp_g_i(g,1);   p2 = 2:(sp_g_i(g,2)-1);
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))'];
            W = R_g_c_p1_p2{g,c}{p1,2}*rxx_m;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ 1 sp_g_i(g,2)-2 sn ]);
%   bound = 3   bound = 4   g = 3
    p1 = sp_g_i(g,1);   p2 = 1;
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))];
            W = R_g_c_p1_p2{g,c}{p1,p2}*rxx_m;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = W;
%   bound = 4   g = 3
    p1 = 2:(sp_g_i(g,1)-1); p2 = 1;
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            rxx_m = [squeeze(unext_g_c_p1_p2_n{gint,c}(p1,p2,:))'];
            W = R_g_c_p1_p2{g,c}{2,p2}*rxx_m;
            unext_g_c_p1_p2_n{gint,c}(p1,p2,:) = reshape(W',[ sp_g_i(g,1)-2 1 sn ]);

end
for g = 1:sg
    unext_g_p1_p2_n{g} = zeros(sp_g_i(g,1), sp_g_i(g,2), sn);
    for c = 1:sc
        unext_g_p1_p2_n{g} = unext_g_p1_p2_n{g} + (1/sc)*unext_g_c_p1_p2_n{g,c};
    end
end

        uprev_g_p1_p2_n = unext_g_c_p1_p2_n;
    end
    for g = 1:sg
        u_g_p1_p2_n_jt{g}(:,:,:,jt) = unext_g_p1_p2_n{g};
    end
end

for g = 1:sg
    K_g_1{g} = zeros(sp_g_i(g,1),sp_g_i(g,2),sjt);
    K_g_2{g} = zeros(sp_g_i(g,1),sp_g_i(g,2),sjt);
end

for g = 1:sg
    K_g_1{g} = -(squeeze(u_g_p1_p2_n_jt{g}(:,:,1,:)).*squeeze(u_g_p1_p2_n_jt{g}(:,:,4,:)) + squeeze(u_g_p1_p2_n_jt{g}(:,:,3,:)).*squeeze(u_g_p1_p2_n_jt{g}(:,:,5,:)));
    K_g_2{g} = -(squeeze(u_g_p1_p2_n_jt{g}(:,:,3,:)).*squeeze(u_g_p1_p2_n_jt{g}(:,:,4,:)) + squeeze(u_g_p1_p2_n_jt{g}(:,:,2,:)).*squeeze(u_g_p1_p2_n_jt{g}(:,:,5,:)));
end

toc;

%%  Графический вывод результатов
x_p1 = cat(2,x_p{1,1},x_p{2,1},x_p{3,1});   
x_p2 = x_p{1,2};
K_p1_p2_jt_1 = cat(1,K_g_1{1},K_g_1{2},K_g_1{3});
K_p1_p2_jt_2 = cat(1,K_g_2{1},K_g_2{2},K_g_2{3});
[sp1,sp2,sj] = size(K_p1_p2_jt_1);

x1min = min(x_p1);                       x1max = max(x_p1);
x2min = min(x_p2);                       x2max = max(x_p2);

Kmin = min(min(min(K_p1_p2_jt_1, [], 1), [], 2), [], 3);   Kmax = max(max(max(K_p1_p2_jt_1, [], 1), [], 2), [], 3);

[X2,X1] = meshgrid(x_p2,x_p1);

[per, sit] = size(T_jt);    clear per;


        for it=1:1:sit
            ZVEC = zeros(3, sp1, sp2);
            ZVEC(1, :, :) = K_p1_p2_jt_1(:,:,it);
            ZVEC(2, :, :) = K_p1_p2_jt_2(:,:,it);
            ZVEC(3, :, :) = 0;
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
            R(1, :, :) = X1;
            R(2, :, :) = X2;
            R(3, :, :) = 0;
            w = typecast(swapbytes(single(R(:))), 'uint8');
            fwrite(f, w);
            fprintf(f, 'CELL_DATA %d\n', (sp1-1) * (sp2-1));
            % No cell data
            fprintf(f, 'POINT_DATA %d\n', sp1 * sp2);
            fprintf(f, 'VECTORS z float\n');
            w = typecast(swapbytes(single(ZVEC(:))), 'uint8');
            fwrite(f, w);
            fclose(f);
            %toc
        end



%%  Формируем матрицы СЛАУ в точках границы 

function [R_g_c_p1_p2,Rd_g_c_p1_p2] = SLAE_boundary()

R_g_c_p1_p2 = cell(sg,sc); Rd_g_c_p1_p2 = cell(sg,sc);
for g = 1:sg
    for c = 1:sc
        R_g_c_p1_p2(g,c) = {cell(sp_g_i(g,1),sp_g_i(g,2))};
        Rd_g_c_p1_p2(g,c) = {cell(sp_g_i(g,1),sp_g_i(g,2))};
    end
end

for c = 1:sc
g = 1;
%   bound = 1   bound = 4   g = 1
    for p1 = 1
        for p2 = 1
            %   Прозрачные граничные условия + Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,1,2};
            WB = [  B_g_c_bound_ind{g,c,4,1};...
                    B_g_c_bound_ind{g,c,1,1}];  
            d_g_p1_p2{g}{p1,p2} = [ zeros(size(B_g_c_bound_ind{g,c,4,1},1),1);...
                                    zeros(size(B_g_c_bound_ind{g,c,1,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 1   g = 1
    for p1 = 1
        for p2 = 2
            %   Прозрачные граничные условия 
            gint = g;   gext = B_g_c_bound_ind{g,c,1,2};
            WB = [B_g_c_bound_ind{g,c,1,1}];     d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{g,c,1,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 1   bound = 2   g = 1
    for p1 = 1
        for p2 = sp_g_i(g,2)
            %   Прозрачные граничные условия + Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            WB = [  B_g_c_bound_ind{g,c,1,1};...
                    B_g_c_bound_ind{g,c,2,1}];  
            d_g_p1_p2{g}{p1,p2} = [ zeros(size(B_g_c_bound_ind{g,c,1,1},1),1);...
                                    zeros(size(B_g_c_bound_ind{g,c,2,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 2   g = 1
    for p1 = 2
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            WB = [B_g_c_bound_ind{gint,c,2,1}]; d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{g,c,2,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 2   bound = 3   g = 1
    for p1 = sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы + Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [  B_g_c_bound_ind{gint,c,2,1}                   zeros(size(B_g_c_bound_ind{gint,c,2,1},1),sn);...
                    zeros(size(B_g_c_bound_ind{gext,c,2,1},1),sn) B_g_c_bound_ind{gext,c,2,1}                  ;...
                    B_g_c_bound_ind{gint,c,3,1}                                                            ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,2,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gext,c,2,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gint,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   g = 1
    for p1 = sp_g_i(g,1)
        for p2 = 2
            %   Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [B_g_c_bound_ind{gint,c,3,1}                                                            ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   bound = 4   g = 1
    for p1 = sp_g_i(g,1)
        for p2 = 1
            %   Условия проскальзывания/слипания + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [  B_g_c_bound_ind{gint,c,4,1},                    zeros(size(B_g_c_bound_ind{gint,c,4,1},1),sn);...
                    zeros(size(B_g_c_bound_ind{gext,c,4,1},1),sn),  B_g_c_bound_ind{gext,c,4,1}                  ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,4,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gext,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 4   g = 1
    for p1 = 2
        for p2 = 1
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            WB = [B_g_c_bound_ind{gint,c,4,1}]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
       
g = 2;
%   bound = 2   g = 2
    for p1 = 2
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            WB = [B_g_c_bound_ind{gint,c,2,1}]; d_g_p1_p2{g}{p1,p2} = [0; 1];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 2   bound = 3   g = 2
    for p1 = sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы + Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [B_g_c_bound_ind{gint,c,2,1}                   zeros(size(B_g_c_bound_ind{gint,c,2,1},1),sn);...
                 zeros(size(B_g_c_bound_ind{gext,c,2,1},1),sn)  B_g_c_bound_ind{gext,c,2,1}                  ;...
                 B_g_c_bound_ind{gint,c,3,1}                                                             ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,2,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gext,c,2,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gint,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   g = 2
    for p1 = sp_g_i(g,1)
        for p2 = 2
            %   Условия проскальзывания/слипания
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [B_g_c_bound_ind{gint,c,3,1}                                                            ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   bound = 4   g = 2
    for p1 = sp_g_i(g,1)
        for p2 = 1
            %   Условия проскальзывания/слипания + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [  B_g_c_bound_ind{gint,c,4,1},                    zeros(size(B_g_c_bound_ind{gint,c,4,1},1),sn);...
                    zeros(size(B_g_c_bound_ind{gext,c,4,1},1),sn),  B_g_c_bound_ind{gext,c,4,1}                  ]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,4,1},1),1);...
                                   zeros(size(B_g_c_bound_ind{gext,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2} zeros(sn,sn)                 ;...
                  zeros(sn,sn)                  C_g_c_p1_p2_m_n{gext,c}{1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 4   g = 2
    for p1 = 2
        for p2 = 1
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            WB = [B_g_c_bound_ind{gint,c,4,1}]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
        
g = 3;
%   bound = 2   g = 3
    for p1 = 2
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            WB = [B_g_c_bound_ind{gint,c,2,1}]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{gint,c,2,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 2   bound = 3   g = 3
    for p1 = sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
            %   Условия свободной границы + Прозрачные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,2,2};
            WB = [  B_g_c_bound_ind{g,c,2,1};...
                    B_g_c_bound_ind{g,c,3,1}];  
            d_g_p1_p2{g}{p1,p2} = [ zeros(size(B_g_c_bound_ind{g,c,2,1},1),1);...
                                    zeros(size(B_g_c_bound_ind{g,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   g = 3
    for p1 = sp_g_i(g,1)
        for p2 = 2
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [B_g_c_bound_ind{g,c,3,1}]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{g,c,3,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 3   bound = 4   g = 3
    for p1 = sp_g_i(g,1)
        for p2 = 1
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,3,2};
            WB = [  B_g_c_bound_ind{g,c,3,1};...
                    B_g_c_bound_ind{g,c,4,1}];     
            d_g_p1_p2{g}{p1,p2} = [ zeros(size(B_g_c_bound_ind{g,c,3,1},1),1);...
                                    zeros(size(B_g_c_bound_ind{g,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
%   bound = 4   g = 3
    for p1 = 2
        for p2 = 1
            %   Прозрачные граничные условия
            gint = g;   gext = B_g_c_bound_ind{g,c,4,2};
            WB = [B_g_c_bound_ind{g,c,4,1}]; 
            d_g_p1_p2{g}{p1,p2} = [zeros(size(B_g_c_bound_ind{g,c,4,1},1),1)];
            WC = [C_g_c_p1_p2_m_n{gint,c}{p1, p2}];
            [R_g_c_p1_p2{g,c}{p1,p2},Rd_g_c_p1_p2{g,c}{p1,p2}] = boundary_solver();
        end
    end
end
function [R,Rd] = boundary_solver()
            [sL,sU,sP,sQ] = lu(sparse(WB)); 
            LB = full(sL);   UB = full(sU);   PB = full(sP);   QB = full(sQ);       clear sL sU sP sQ; 
            PB1 = inv(LB(1:rank(LB),:))*PB(1:rank(LB),:);
            
            [R,jrow] = rref(UB');   PU_k = zeros(length(jrow),size(UB,1),size(UB,1));   PU = eye(size(UB,1));
            for k = 1:length(jrow)
                PU_k(k,:,:) = eye(size(UB,1));
                PU_k(k,k,:) = 0;        PU_k(k,jrow(k),:) = 0;
                PU_k(k,k,jrow(k)) = 1;  PU_k(k,jrow(k),k) = 1;
                PU = squeeze(PU_k(k,:,:))*PU;
            end
            clear PU_k;
            [R,jcol] = rref(UB);   QU_k = zeros(length(jcol),size(UB,2),size(UB,2));   QU = eye(size(UB,2));
            for k = 1:length(jcol)
                QU_k(k,:,:) = eye(size(UB,2));
                QU_k(k,:,k) = 0;        QU_k(k,:,jcol(k)) = 0;
                QU_k(k,jcol(k),k) = 1;  QU_k(k,k,jcol(k)) = 1;
                QU = QU*squeeze(QU_k(k,:,:));
            end
            clear QU_k;
            WUB = PU*UB*QU; WPB1 = PU*PB1;  clear PU UB PB1;
            UB1 = WUB(1:rank(WUB),1:rank(WUB));   UB2 = WUB(1:rank(WUB),rank(WUB)+1:end);
            WW2 = inv(UB1)*UB2; 
            PB2 = inv(UB1)*WPB1(1:rank(WUB),:);   clear UB1 UB2;
            WW = QB*QU*[-WW2; eye(size(WUB,2)-rank(WUB))];   clear WW2;
            Wd = QB*QU*[PB2; zeros(size(WUB,2)-rank(WUB),size(WB,1))];    clear PB2;
            
%             [sL,sU,sP,sQ] = lu(sparse(WC*WW)); 
%             LC = full(sL);   UC = full(sU);   PC = full(sP);   QC = full(sQ);       clear sL sU sP sQ;
%             PC1 = inv(LC(1:rank(LC),:))*PC(1:rank(LC),:);
%             
%             [R,jrow] = rref(UC');   PU_k = zeros(length(jrow),size(UC,1),size(UC,1));   PU = eye(size(UC,1));
%             for k = 1:length(jrow)
%                 PU_k(k,:,:) = eye(size(UC,1));
%                 PU_k(k,k,:) = 0;        PU_k(k,jrow(k),:) = 0;
%                 PU_k(k,k,jrow(k)) = 1;  PU_k(k,jrow(k),k) = 1;
%                 PU = squeeze(PU_k(k,:,:))*PU;
%             end
%             clear PU_k;
%             [R,jcol] = rref(UC);   QU_k = zeros(length(jcol),size(UC,2),size(UC,2));   QU = eye(size(UC,2));
%             for k = 1:length(jcol)
%                 QU_k(k,:,:) = eye(size(UC,2));
%                 QU_k(k,:,k) = 0;        QU_k(k,:,jcol(k)) = 0;
%                 QU_k(k,jcol(k),k) = 1;  QU_k(k,k,jcol(k)) = 1;
%                 QU = QU*squeeze(QU_k(k,:,:));
%             end
%             clear QU_k;
%             
%             WUC = PU*UC*QU; WPC1 = PU*PC1;
%             R =WW*QC*QU*inv(WUC(1:rank(WUC),:))*WPC1(1:rank(WUC),:);
%             Rd = (eye(size(R,1))-R*WC)*Wd;


            [LC,UC,PC] = lu(WC*WW);
            R = WW*inv(UC)*inv(LC(1:rank(LC),:))*PC(1:rank(LC),:);
            Rd = (eye(size(R,1))-R*WC)*Wd;
end
end

%% Матрицы соответствующие внутреним характеристикам

function C_g_c_p1_p2_m_n = internal_characteristic()
    
C_g_c_p1_p2_m_n = cell(sg,sc);
for g = 1:sg
    for c = 1:sc
        C_g_c_p1_p2_m_n(g,c) = {cell(sp_g_i(g,1),sp_g_i(g,2))};
    end
end


for g = 1:sg
    for c = 1:sc

%         for p1 = 2:(sp_g_i(g,1)-1)
%             for p2 = 2:(sp_g_i(g,2)-1)
%                 C_g_c_p1_p2_m_n{g,c}(p1,p2) = {eye(sn)};
%             end
%         end

    for p1 = 1
        for p2 = 1
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)<=tollambd),find(lambd_g_i_k(g,2,:)<=tollambd),:,:),1),2))};
        end
    end
    
    for p1 = 1
        for p2 = 2
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)<=tollambd),:,:,:),1),2))};
        end
    end
    
    for p1 = 1
        for p2 = sp_g_i(g,2)
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)<=tollambd),find(lambd_g_i_k(g,2,:)>=-tollambd),:,:),1),2))};
        end
    end
    
    for p1 = 2
        for p2 = sp_g_i(g,2)
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(:,find(lambd_g_i_k(g,2,:)>=-tollambd),:,:),1),2))};
        end
    end
    
    for p1 = sp_g_i(g,1)
        for p2 = sp_g_i(g,2)
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)>=-tollambd),find(lambd_g_i_k(g,2,:)>=-tollambd),:,:),1),2))};
        end
    end
    
    for p1 = sp_g_i(g,1)
        for p2 = 2
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)>=-tollambd),:,:,:),1),2))};
        end
    end
    
    for p1 = sp_g_i(g,1)
        for p2 = 1
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)>=-tollambd),find(lambd_g_i_k(g,2,:)<=tollambd),:,:),1),2))};
        end
    end
    
    for p1 = 2
        for p2 = 1
C_g_c_p1_p2_m_n{g,c}(p1,p2) = {squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(:,find(lambd_g_i_k(g,2,:)<=tollambd),:,:),1),2))};
        end
    end
end
end
end

%%  Условия сопряжения на границах

function COND_g_c_bound_ind = boundary_conditions()
    
%   Матрица свободных граничных условий
Bfree_i_m_n = zeros(sdim,2,sn);    
Bfree_i_m_n(1,1,3) = 1; Bfree_i_m_n(1,2,1) = 1;
Bfree_i_m_n(2,1,3) = 1; Bfree_i_m_n(2,2,2) = 1;

%   Матрица  условий слипания
Badhes_i_m_n = zeros(sdim,4,2*sn); 
Badhes_i_m_n(1,1,1) = 1;    Badhes_i_m_n(1,1,1+sn) = -1;
Badhes_i_m_n(1,2,3) = 1;    Badhes_i_m_n(1,2,3+sn) = -1;
Badhes_i_m_n(1,3,4) = 1;    Badhes_i_m_n(1,3,4+sn) = -1;
Badhes_i_m_n(1,4,5) = 1;    Badhes_i_m_n(1,4,5+sn) = -1;

Badhes_i_m_n(2,1,2) = 1;    Badhes_i_m_n(2,1,2+sn) = -1;
Badhes_i_m_n(2,2,3) = 1;    Badhes_i_m_n(2,2,3+sn) = -1;
Badhes_i_m_n(2,3,4) = 1;    Badhes_i_m_n(2,3,4+sn) = -1;
Badhes_i_m_n(2,4,5) = 1;    Badhes_i_m_n(2,4,5+sn) = -1;

%   Матрица  условий проскальзывания
Bsliding_i_m_n = zeros(sdim,4,2*sn); 
Bsliding_i_m_n(1,1,1)    = 1;   Bsliding_i_m_n(1,1,1+sn) = -1;
Bsliding_i_m_n(1,2,3)    = 1;    
Bsliding_i_m_n(1,3,3+sn) = 1;    
Bsliding_i_m_n(1,4,4)    = 1;   Bsliding_i_m_n(1,4,4+sn) = -1;
% Bsliding_i_m_n(1,5,2)    = 1;   Bsliding_i_m_n(1,5,2+sn) = -1;


Bsliding_i_m_n(2,1,2)    = 1;   Bsliding_i_m_n(2,1,2+sn) = -1;
Bsliding_i_m_n(2,2,3)    = 1;    
Bsliding_i_m_n(2,3,3+sn) = 1;    
Bsliding_i_m_n(2,4,5)    = 1;   Bsliding_i_m_n(2,4,5+sn) = -1;
% Bsliding_i_m_n(2,5,1)    = 1;   Bsliding_i_m_n(2,5,1+sn) = -1;

COND_g_c_bound_ind = cell(sg,sc,4,2);
for c = 1:sc
g = 1;  
    bound = 1;
    %   Прозрачные граничные условия на границе p1 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)>tollambd),:,:,:),1),2));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
    bound = 2;
    %   Свободные граничные условия на границе p2 = sp_g_i(g,2)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Bfree_i_m_n(2,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
    bound = 3;
    %   Условия проскальзывания/слипания на границе p1 = sp_g_i(g,1)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Badhes_i_m_n(1,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 2;
    bound = 4;
    %   Прозрачные граничные условия на границе p2 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(:,find(lambd_g_i_k(g,2,:)>tollambd),:,:),1),2));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
g = 2;  
    bound = 1;
    %   Условия проскальзывания/слипания на границе p1 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Badhes_i_m_n(1,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 1;
    bound = 2;
    %   Свободные граничные условия на границе p2 = sp_g_i(g,2)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Bfree_i_m_n(2,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
    bound = 3;
    %   Условия проскальзывания/слипания на границе p1 = sp_g_i(g,1)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Badhes_i_m_n(1,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 3;
    bound = 4;
    %   Прозрачные граничные условия на границе p2 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(:,find(lambd_g_i_k(g,2,:)>tollambd),:,:),1),2));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
g = 3;  
    bound = 1;
    %   Условия проскальзывания/слипания на границе p1 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Badhes_i_m_n(1,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 2;
    bound = 2;
    %   Свободные граничные условия на границе p2 = sp_g_i(g,2)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(Bfree_i_m_n(2,:,:));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
    bound = 3;
    %   Прозрачные граничные условия на границе p1 = sp_g_i(g,1)
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(find(lambd_g_i_k(g,1,:)<-tollambd),:,:,:),1),2));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
    bound = 4;
    %   Прозрачные граничные условия на границе p2 = 1
    COND_g_c_bound_ind{g,c,bound,1} = squeeze(sum(sum(C_g_c_k1_k2_m_n{g,c}(:,find(lambd_g_i_k(g,2,:)>tollambd),:,:),1),2));
    COND_g_c_bound_ind{g,c,bound,2} = 0;
end
end


%%  Формируем базисные полиномы
function H_g_i_pw_k_p = Basic_polynoms()
    
%%   Вспомагательные процедуры
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

H_g_i_pw_k_p = cell(sg,sdim);
for g = 1:sg
    for i = 1:sdim
        H_g_i_pw_k_p{g,i} = zeros(sp_g_i(g,i),sn,sp_g_i(g,i));
    end
end

for g = 1:sg
    for i = 1:sdim
        for pw = 1
            for k = find(lambd_g_i_k(g,i,:)<=tollambd)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(1,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_g_i_k(g,i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw+1) = polyval(pbase_l_k(2,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_g_i_k(g,i,k)*deltat));
            end
        end

        for pw = 2:(sp_g_i(g,i)-1)
            for k = find(lambd_g_i_k(g,i,:)<=tollambd)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(1,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_g_i_k(g,i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw+1) = polyval(pbase_l_k(2,:), (-1-(2/(h_g_i_p{g,i}(pw)))*lambd_g_i_k(g,i,k)*deltat));
            end
            for k = find(lambd_g_i_k(g,i,:)>=-tollambd)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(2,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_g_i_k(g,i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw-1) = polyval(pbase_l_k(1,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_g_i_k(g,i,k)*deltat));
            end
        end

        for pw = sp_g_i(g,i)
            for k = find(lambd_g_i_k(g,i,:)>=-tollambd)
H_g_i_pw_k_p{g,i}(pw,k,pw)   = polyval(pbase_l_k(2,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_g_i_k(g,i,k)*deltat));
H_g_i_pw_k_p{g,i}(pw,k,pw-1) = polyval(pbase_l_k(1,:), (1-(2/(h_g_i_p{g,i}(pw-1)))*lambd_g_i_k(g,i,k)*deltat));
            end
        end
    end
end

for g = 1:sg
    for i = 1:sdim
        H_g_i_pwk_p{g,i} = sparse(reshape(H_g_i_pw_k_p{g,i}, [ sp_g_i(g,i)*sn sp_g_i(g,i)]));
    end
end
end

end