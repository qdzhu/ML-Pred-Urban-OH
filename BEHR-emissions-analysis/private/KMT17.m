function val = KMT17(TEMP, M)
K170 = 5.0D-30*M*(TEMP/300).^-1.5;
K17I = 1.0D-12;
KR17 = K170 ./ K17I;
FC17 = 0.17*exp(-51/TEMP)+exp(-TEMP/204);
NC17 = 0.75-1.27*(log10(FC17));
F17 = 10.^(log10(FC17)/(1.0+(log10(KR17)/NC17).^2));
val = (K170*K17I*F17)/(K170+K17I);
end