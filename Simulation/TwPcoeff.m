function coeff = TwPcoeff(n,cori)
%short function to spit out the coefficient to multiply TwP by in the
%growth sim to minimize twinning in the c-axis growth direction.

if n == 1 || n == 5
    coeff = cosd(cori)^2;
    return
end

if n == 2 || n == 6
    coeff = sind(cori-135)^2;
    return
end

if n == 3 || n == 7
    coeff = sind(cori)^2;
    return
end

if n == 4 || n == 8
    coeff = sind(cori-45)^2;
    return
end