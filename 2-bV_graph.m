function B = bV_graph(V,l,a,b,n1,n2,n3)
%BV_GRAPH функциZ построения графика зависимости нормализованной частоты B 
% от нормализованной частоты V для диэлектрического волновода W-типа
% входные данные:
% V  -  вектор значений нормализованной частоты
% a   -  внутренний радиус внутренней оболочки
% b  -  внешний радиус внутренней оболочки
% n1 - показатель преломления волноведущей жилы
% n2 - показатель преломления средней оболочки
% n3 - показатель преломления внешней оболочки

    e1 = n1.^2;
    e2 = n2.^2;
    e3 = n3.^2;
    
    nn = 1:1000;
    B  = 1./(exp((nn-500)./25)+1); % сетка предварительных значений
    
    lB = length(B);
    lV = length(V);
    
    F = zeros(lB, lV);
    
    for ii = 1:lB
        for jj = 1:lV
            k = V(jj)./(a*sqrt(e1-e3)); % это омега деленная на скорость света
            betta = sqrt((n1.^2-n3.^2)*B(ii)+n3.^2).*k;
            u = sqrt(e1.*k.^2-betta.^2);
            w = sqrt(betta.^2-e2.*k.^2);
            v = sqrt(betta.^2-e3.*k.^2);

            Jlua = besselj(l,u.*a);
            Klwa = besselk(l,w.*a);
            Ilwa = besseli(l,w.*a);
            Klwb = besselk(l,w.*b);
            Ilwb = besseli(l,w.*b);
            Klvb = besselk(l,v.*b);

            J_lua = -besselj(l+1,u.*a)+l.*Jlua./(u.*a);
            K_lwa = -besselk(l+1,w.*a)+l.*Klwa./(w.*a);
            I_lwa = -besseli(l+1,w.*a)+l.*Ilwa./(w.*a);
            K_lwb = -besselk(l+1,w.*b)+l.*Klwb./(w.*b);
            I_lwb = -besseli(l+1,w.*b)+l.*Ilwb./(w.*b);
            K_lvb = -besselk(l+1,v.*b)+l.*Klvb./(v.*b);
            
            
            buw = betta.*l.*(u.^(-2)+w.^(-2))./a;
            bwv = betta.*l.*(w.^(-2)-v.^(-2))./b;
            ke1 = i.*e1.*k./u;
            ke2 = i.*e2.*k./w;
            ke3 = i.*e3.*k./v;
            
            kw  = i.*k./w;
            kv  = i.*k./v;
            ku  = i.*k./u;
            
    
            S = [Jlua,        -Klwa,     - Ilwa,         0,          0,          0,       0,         0;
                   0,         -Klwb,      -Ilwb,        Klvb,        0,          0,       0,         0;
                -ke1.*J_lua, -ke2.*K_lwa, -ke2.*I_lwa,   0,        buw.*Jlua,    0,       0,         0;
                   0,        -ke2.*K_lwb, -ke2.*I_lwb, ke3.*K_lvb,   0,          0,       0,       bwv.*Klvb;
                buw.*Jlua,      0,           0,          0,       ku.*J_lua, kw.*K_lwa, kw.*I_lwa,   0;
                   0,           0,           0,         bwv.*Klvb,    0,     kw.*K_lwb, kw.*I_lwb, -kv.*K_lvb;
                   0,           0,           0,          0,           0,       -Klwb,    -Ilwb,     Klvb;
                   0,           0,           0,          0,         Jlua,      -Klwa,    -Ilwa,      0];
            
            F(ii, jj) = (det(S));
        end
    end
    figure;
    hold on;
    grid on;
    contour(V, B, F,[0 0]);
    
end
