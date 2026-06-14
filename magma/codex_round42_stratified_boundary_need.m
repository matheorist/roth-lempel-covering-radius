// codex_round42_stratified_boundary_need.m
// Counts boundary classes after removing those certified theoretically by
// q-21+3r > 0, then by an explicit unshifted curve witness.
SetEchoInput(false);
function K(x) return Sprint(x); end function;
function CK(S,d) return Join(Sort([K(x):x in S]),",") cat "|" cat K(d); end function;
function CR(F,S,d)
 first:=true;k0:="";S0:=S;d0:=d;
 for x in S do for y in S do if x ne y then
  a:=(y-x)^-1;b:=-a*x;T:={a*z+b:z in S};e:=a*d+3*b;k:=CK(T,e);
  if first or k lt k0 then first:=false;k0:=k;S0:=SetToSequence(T);d0:=e;end if;
 end if;end for;end for;return k0,S0,d0;
end function;
function Reps(F)
 A:={x:x in F|x ne 0 and x ne 1}; D:=AssociativeArray(); R:=[];
 for X in Subsets(A,4) do for d in F do k,S,e:=CR(F,{F!0,F!1} join X,d);
  if not IsDefined(D,k) then D[k]:=true;Append(~R,<S,e>);end if; end for;end for;return R;
end function;
function Rv(S,d) return #[I:I in Subsets({1..6},3)|&+[S[i]:i in SetToSequence(I)] eq d];end function;
function Bset(S,d)
 th:=&+S-d; B:=SequenceToSet(S);
 for I in Subsets({1..6},2) do J:=SetToSequence(I);Include(~B,th-S[J[1]]-S[J[2]]);end for;return B;
end function;
procedure Scan(q)
 F<w>:=GF(q);total:=0; theorem:=0; exactCurveOnly:=0; residual:=0; residualByR:=AssociativeArray();
 for X in Reps(F) do
  total+:=1; r:=Rv(X[1],X[2]);
  if q-21+3*r gt 0 then theorem+:=1;
  elif #Bset(X[1],X[2]) lt q then exactCurveOnly+:=1;
  else
   residual+:=1;k:=Sprint(r);if not IsDefined(residualByR,k) then residualByR[k]:=0;end if;residualByR[k]+:=1;
  end if;
 end for;
 printf "STRATIFIED_BOUNDARY q=%o total=%o theoremByR=%o extraCurve=%o residualCertificates=%o residualByR=%o\n",
  q,total,theorem,exactCurveOnly,residual,residualByR;
end procedure;
for q in [13,16,17,19] do Scan(q);end for;quit;
