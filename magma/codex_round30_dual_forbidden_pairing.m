// codex_round30_dual_forbidden_pairing.m
// Compares twisted-cubic forbidden-set sizes for delta and theta=sum(S)-delta.
SetEchoInput(false);
function K(x) return Sprint(x);end function;
function CK(S,d) return Join(Sort([K(x):x in S]),",") cat "|" cat K(d);end function;
function CR(F,S,d)
 first:=true;k0:="";S0:=S;d0:=d;
 for x in S do for y in S do if x ne y then
  a:=(y-x)^-1;b:=-a*x;T:={a*z+b:z in S};e:=a*d+3*b;k:=CK(T,e);
  if first or k lt k0 then first:=false;k0:=k;S0:=SetToSequence(T);d0:=e;end if;
 end if;end for;end for;return k0,S0,d0;
end function;
function R(F)
 A:={x:x in F|x ne 0 and x ne 1};D:=AssociativeArray();Z:=[];
 for X in Subsets(A,4) do for d in F do k,S,e:=CR(F,{F!0,F!1} join X,d);
  if not IsDefined(D,k) then D[k]:=true;Append(~Z,<S,e>);end if;end for;end for;return Z;
end function;
function Bs(S,d)
 th:=&+S-d;B:=SequenceToSet(S);
 for I in Subsets({1..6},2) do J:=SetToSequence(I);Include(~B,th-S[J[1]]-S[J[2]]);end for;return #B;
end function;
procedure Scan(q)
 F<w>:=GF(q);sameCount:=0;differentCount:=0;maxgap:=0;
 for X in R(F) do S:=X[1];d:=X[2];th:=&+S-d;a:=Bs(S,d);b:=Bs(S,th);
  if a eq b then sameCount+:=1; else differentCount+:=1;if Abs(a-b) gt maxgap then maxgap:=Abs(a-b);end if;end if;
 end for;
 printf "DUAL_B q=%o same=%o different=%o maxGap=%o\n",q,sameCount,differentCount,maxgap;
end procedure;
for q in [7,8,9,11,13,16,17,19,23] do Scan(q);end for;quit;
