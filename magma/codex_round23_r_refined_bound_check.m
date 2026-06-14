// codex_round23_r_refined_bound_check.m
// Checks |B_{S,delta}| <= 21-3 r_S(delta) and records equality cases.

SetEchoInput(false);
function EK(x) return Sprint(x); end function;
function CK(S,d) return Join(Sort([EK(x):x in S]),",") cat "|" cat EK(d); end function;
function CR(F,S,d)
    first:=true; k0:=""; T0:=S; d0:=d;
    for x in S do for y in S do if x ne y then
        a:=(y-x)^-1; b:=-a*x; T:={a*z+b:z in S}; e:=a*d+3*b; k:=CK(T,e);
        if first or k lt k0 then first:=false; k0:=k; T0:=SetToSequence(T); d0:=e; end if;
    end if; end for; end for;
    return k0,T0,d0;
end function;
function Reps(F)
    A:={x:x in F|x ne F!0 and x ne F!1}; D:=AssociativeArray(); R:=[];
    for X in Subsets(A,4) do for d in F do
        k,S,e:=CR(F,{F!0,F!1} join X,d);
        if not IsDefined(D,k) then D[k]:=true; Append(~R,<S,e>); end if;
    end for; end for; return R;
end function;
function RVal(S,d) return #[I:I in Subsets({1..6},3)|&+[S[i]:i in SetToSequence(I)] eq d]; end function;
function BSize(S,d)
    th:=&+S-d; B:=SequenceToSet(S);
    for I in Subsets({1..6},2) do J:=SetToSequence(I); Include(~B,th-S[J[1]]-S[J[2]]); end for;
    return #B;
end function;
procedure Scan(q)
    F<w>:=GF(q); failures:=0; equal:=0; maxByR:=AssociativeArray();
    for X in Reps(F) do
        r:=RVal(X[1],X[2]); b:=BSize(X[1],X[2]); bd:=21-3*r;
        if b gt bd then failures+:=1; end if;
        if b eq bd then equal+:=1; end if;
        k:=Sprint(r); if not IsDefined(maxByR,k) or b gt maxByR[k] then maxByR[k]:=b; end if;
    end for;
    printf "R_BOUND q=%o failures=%o equalityClasses=%o maxima=%o\n",q,failures,equal,maxByR;
end procedure;
for q in [7,8,9,11,13,16,17,19,23] do Scan(q); end for;
quit;
