function [Ex,Ey,Ez,Hx,Hy,Hz,beta] = my_field_distrib(lmbd,sst,l,m,a,b,n1,n2,n3,x,y)  
%Работает для r <= a%
%phi отсчитывается от вертикали по часовой%

r = sqrt(x.^2 + y.^2);
dn = (n1-n3).*(n1+n3);
wc = 2.*pi./lmbd;
V = wc.*a.*sqrt(dn);
c = b./a;
Bb = bV_dispersion_law(V,l,m,sst,c,n1,n2,n3);
bn = Bb(end);

N = sqrt(bn.*dn+n3.^2);
beta = N.*wc;

if bn ~= 0
    u = V.*sqrt(1-bn)./a;
    w = wc.*sqrt(bn.*dn + n3^2 - n2^2);
    v = V.*sqrt(bn)./a;
else
    disp('The normalize frequency V is lower then cutoff frequency. ');
    Ex = 0; Ey = 0; Ez = 0;
    Hx = 0; Hy = 0; Hz = 0;
    return;
end

Jua = besselj(l, u.*a);
Jwa = besselj(l, w.*a);
Jur = besselj(l, u.*r);
Kwa = besselk(l, w.*a);
Kwb = besselk(l, w.*b);
Kvb = besselk(l, v.*b);
Kwr = besselk(l, w.*r);
Kvr = besselk(l, v.*r);
Iwa = besseli(l, w.*a);
Iwb = besseli(l, w.*b);
Iwr = besseli(l, w.*r);
bl = beta.*l;
iwc = 1i.*wc;
Jua1 = -besselj(l + 1, u.*a) + l.*Jua./(u.*a);  
Jwa1 = -besselj(l + 1, w.*a) + l.*Jwa./(w.*a);
Kwa1 = -besselk(l + 1, w.*a) + l.*Kwa./(w.*a);
Kwb1 = -besselk(l + 1, w.*b) + l.*Kwb./(w.*b);
Kvb1 = -besselk(l + 1, v.*b) + l.*Kvb./(v.*b);
Iwa1 = -besseli(l + 1, w.*a) + l.*Iwa./(w.*a);
Iwb1 = -besseli(l + 1, w.*b) + l.*Iwb./(w.*b);
Jur1 = -besselj(l + 1, u.*r) + l.*Jur./(u.*r); 
Kwr1 = -besselk(l + 1, w.*r) + l.*Kwr./(w.*r);
Iwr1 = -besseli(l + 1, w.*r) + l.*Iwr./(w.*r);
Kvr1 = -besselk(l + 1, v.*r) + l.*Kvr./(v.*r);


S = [Jua,    -Kwa,    -Iwa,    0,    0,    0,    0,    0;
     0,    Kwb,    Iwb,    Kvb,    0,    0,    0,    0;
     0,    0,    0,    0,    Jua,    -Kwa,    -Iwa,    0;
     0,    0,    0,    0,    0,    Kwb,    Iwb,    -Kvb;
     bl./(u.^2.*a).*Jua,    bl./(w.^2.*a).*Kwa,    bl./(w.^2.*a).*Iwa,    0,    iwc./u.*Jua1,    iwc./w.*Kwa1,    iwc./w.*Iwa1,    0;
     0,    bl./(w.^2.*b).*Kwb,    bl./(w.^2.*b).*Iwb,    -bl./(v.^2.*b).*Kvb,    0,    iwc./w.*Kwb1,    iwc./w.*Iwb1,    -iwc./v.*Kvb1;
     -iwc.*n1.^2./u.*Jua1,    -iwc.*n2.^2./w.*Kwa1,    -iwc.*n2.^2./w.*Jwa1,    0,    bl./(u.^2.*a).*Jua,    bl./(w.^2.*a).*Kwa,    bl./(w.^2.*a).*Iwa,    0;
     0,    -iwc.*n2.^2./w.*Kwb1,    -iwc.*n2.^2./w.*Iwb1,    iwc.*n3.^2./v.*Kvb1,    0,    bl./(w.^2.*b).*Kwb,    bl./(w.^2.*b).*Iwb,    -bl./(v.^2.*b).*Kvb];

if l ~= 0
        coeff = null(S);
        factor = 1/coeff(1);
        A1 = 1;
        A2 = factor.*coeff(2);
        A3 = factor.*coeff(3);
        A4 = factor.*coeff(4);
        B1 = factor.*coeff(5);
        B2 = factor.*coeff(6);
        B3 = factor.*coeff(7);
        B4 = factor.*coeff(8);
    
else
    if sst == 'EH'
        A1 = 0;
        A4 = 0;
        A2 = 1;
        A3 = -Klwr./Ilwr.*A1;
        B1 = 1;
    else
        B1 = 0;
        B4 = 0;
        B2 = 1;  %т.к не зависит от А
        B3 = -Klwr./Ilwr.*B1;
        A1 = 1;
    end
end

if r <= a
    Ez = abs(A1.*Jur);
    Ex = abs(1i./(u.^2).*(A1.*beta.*Jur1.*u + iwc.*B1.*l./r.*Jur)).*(x./sqrt(x.^2 + y.^2)) + abs(1i./(u.^2).*(1i.*A1.*beta.*l./r.*Jur -B1.*wc.*Jur1.*u)).*(y./sqrt(x.^2 + y.^2));
    Ey = abs(1i./(u.^2).*(A1.*beta.*Jur1.*u + iwc.*B1.*l./r.*Jur)).*(y./sqrt(x.^2 + y.^2)) + abs(1i./(u.^2).*(1i.*A1.*beta.*l./r.*Jur -B1.*wc.*Jur1.*u)).*(x./sqrt(x.^2 + y.^2));
    Hx = 0;
    Hy = 0;
    Hz = 0;
    return
elseif (r >= a) && (r <= b)
    return
elseif r >= b
    return
end
end
