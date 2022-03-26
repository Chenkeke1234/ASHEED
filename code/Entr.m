function EntranPCH =Entr(a,b)
if b<=87
    EntranPCH=50*10^(-9)*a+10*10^(-12)*a*b^2;
else
    EntranPCH=50*10^(-9)*a+1.3*10^(-15)*a*b^4;
end

end

