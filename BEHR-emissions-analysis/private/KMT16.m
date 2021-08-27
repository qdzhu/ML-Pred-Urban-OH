function val = KMT16(TEMP, M)
K160 = 8D-27*M*(TEMP/300) .^ -3.5;
K16I = 3.0D-11*(TEMP/300) .^ -1;
KR16 = K160/K16I;
FC16 = 0.5;
NC16 = 0.75-1.27*(log10(FC16));
F16 = 10 .^ (log10(FC16)/(1+(log10(KR16)/NC16) .^ 2));
val = (K160*K16I)*F16/(K160+K16I);
end