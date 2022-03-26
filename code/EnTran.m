function EntranPCH=EnTran(a,b,c,d)
D=7569;
if d<=D
   EntranPCH=a*c+10*0.000000000001*c*d;
else
   EntranPCH=a*c+c*0.0013*0.000000000001*d^2;
end
end

