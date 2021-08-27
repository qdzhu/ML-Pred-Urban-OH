function val = KMT15(TEMP, M)
k150 = 8.6D-29*M*(TEMP/300) .^ -3.1;
k15i = 9.0D-12*(TEMP/300) .^ -0.85;
kr15 = k150 ./ k15i;
fc15 = 0.48;
nc15 = 0.75-1.27*(log10(fc15));
f15 = 10 .^ (log10(fc15)/(1+(log10(kr15)/nc15).^2));
val = (k150*k15i)*f15/(k150+k15i);
end