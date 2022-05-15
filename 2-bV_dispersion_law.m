lambda = [0.5:0.01:10];
n1 = 1.445;
n2 = 1.435;
n3 = 1.440;

a = 10;
b = 13;
c = b./a;

sst = 'EH';
m_max = 5;

l = 1;
V = 2.*pi./lambda.*a*sqrt(n1.^2-n3.^2);

Bb = bV_dispersion_law(V, l, m_max, sst, c, n1, n2, n3);

figure
grid on;
hold on;
plot(V,Bb);

function Bb = bV_dispersion_law(Vin, lin, m_max, sst, cin, n1, n2, n3)

global V l e1 e2 e3 c
    l = lin;
    e1 = n1.^2;
    e2 = n2.^2;
    e3 = n3.^2;
    c = cin;
    
    nn = 1:1000;
    b  = 1./(exp((nn-500)./25)+1); % сетка предварительных значений

    lV = length(Vin);
    
    switch sst
    case 'EH'
        fn = @EH; % function handle of chracteristic function of EH modes
    case 'HE'
        fn = @HE; % function handle of chracteristic function  of HE modes
    otherwise
        disp('Unknow type of mode. Only HE and EH modes are possible.');
        Bb = NaN;
        return;
    end
    
    Bb = zeros(m_max,lV);
    bv = zeros(m_max,1);
    for k = 1:lV
        V = Vin(k);
        F  = fn(b);
        F2 = F(1:end-1).*F(2:end); 
        inds  = find(F2<0);
        len_i = length(inds);
        if (len_i>m_max) 
            len_i = m_max;
        end
        for ii = 1:len_i
            bv(ii) = fzero(fn,[b(inds(ii)),b(inds(ii)+1)]);
        end
        Bb(1:len_i,k) = bv(1:len_i); 
    end
    
end

function [E,H] = W_charact_function(B)
%W_CHARACT_FUNCTION функция левой части характеристического уравнения для
% волновода W-типа

global V l e1 e2 e3 c % значения этих параметров должны быть переданы через область глобальной видимости переменных

% V - нормализованная частота волновода
% l - азимутальное число
% e1 - диэлектрическая проницаемость жилы
% e2 - диэлектрическая проницаемость средней оболочки
% e3 - диэлектрическая проницаемость внешней оболочки
% c - отношение внешнего радиуса средней оборочки к внутреннему (радиусу волноведущей жилы)

% далее следуют промежуточные вычисления
	gamma   = (e3-e2)/(e1-e3);
	delta   = e3/(e1-e3);
	ksita   = e2/(e1-e3);

	sqB 	= sqrt(B);
	sq_B	= sqrt(1-B);
	sqBg	= sqrt(B+gamma);

	ua  	= V.*sq_B;
	wa  	= V.*sqBg;
	wb  	= V.*c.*sqBg;
	vb  	= V.*c.*sqB;

	Jlua	= besselj(l,ua);
	J_lua   = besselj(l+1,ua);

	Klwa	= besselk(l,wa);
	K_lwa   = besselk(l+1,wa);  
	Ilwa	= besseli(l,wa);   
	I_lwa   = besseli(l+1,wa);  
	Klwb	= besselk(l,wb);
	K_lwb   = besselk(l+1,wb);  
	Ilwb	= besseli(l,wb);   
	I_lwb   = besseli(l+1,wb);  

	Z   	= besselk(l+1,vb)./besselk(l,vb);

	alfa1   = (K_lwa./Klwb-I_lwa./Ilwb)./(Klwa./Klwb-Ilwa./Ilwb);
	alfa2   = (K_lwb./Klwb-I_lwb./Ilwb)./(Klwa./Klwb-Ilwa./Ilwb);
	beta1   = (I_lwa./Ilwa-K_lwa./Klwa)./(Ilwb./Ilwa-Klwb./Klwa);
	beta2   = (I_lwb./Ilwa-K_lwb./Klwa)./(Ilwb./Ilwa-Klwb./Klwa);

	clear Klwa K_lwa Klwb K_lwb Ilwa I_lwa Ilwb I_lwb;

	p   	= e2.*B.*sqBg.*beta2-e3.*(B+gamma).*sqB.*Z;
	r   	= B.*sqBg.*beta2-(B+gamma).*sqB.*Z;
	h   	= e2.*(1-B).*sqBg.*Jlua.*alfa1+e1.*(B+gamma).*sq_B.*J_lua;
	g   	= (1-B).*sqBg.*Jlua.*alfa1+(B+gamma).*sq_B.*J_lua;

	clear Z;

	K1  	= V.*c.*p.*r+l.*((e3-e2).*(B+delta).*r+gamma.*p);
	JikB1   = Jlua.*alfa2.*beta1.*B.*(1-B).*(B+gamma);

	qb  	= JikB1.*(e1-e3)./K1;
	qc  	= (c.*(V.*JikB1.*e2).^2)./K1;
	d1  	= (l.*gamma+r.*c.*V).*V.*e2.^2./(e1-e3);
	d2  	= (l.*gamma.*(B+delta)+p.*c.*V./(e1-e3)).*V;
	f1  	= Jlua.*l.*((e1-e2).*B+(e1-e2).*delta);
	f2  	= Jlua.*l.*(gamma+1);
	dhg 	= Jlua.*(1-B).*sqBg.*alfa1.*(e1-e2);
	A   	= Jlua.*l.*(gamma+1).*(l.*(e3-e2).*((B+delta+ksita).^2)+V.*c.*((B+delta).*p+r.*ksita.*e2));

	aa  	= V.*e1;
	bb  	= -(dhg.*V + f1 + qb.*d1 + e1.*( f2 + qb.*d2));
	cc  	= qc + qb.*A + dhg.*(f2 + qb.*d2);

	Det 	= bb.^2-4.*aa.*cc;
    
	% вычисление целевых значений
	H   	= g - (-bb + sqrt(Det))./(2.*aa);
	E   	= g - (-bb - sqrt(Det))./(2.*aa);

end

function H = HE(B)
       [a,H] = W_charact_function(B);
 end
 function E = EH(B)
       [E,a] = W_charact_function(B);
 end
