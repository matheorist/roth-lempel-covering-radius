// codex_round16_m2_formula_boundary_codex.m
// Verification that m_2(C)=28(q-1)-3 r_S(delta) on boundary fields.

SetEchoInput(false);

function EK(x) return Sprint(x); end function;
function CK(S,d) return Join(Sort([EK(x):x in S]),",") cat "|" cat EK(d); end function;
function CR(F,S,d)
    first:=true; bk:=""; bS:=S; bd:=d;
    for x in S do for y in S do if x ne y then
        a:=(y-x)^-1; b:=-a*x; T:={a*z+b:z in S}; e:=a*d+3*b; k:=CK(T,e);
        if first or k lt bk then first:=false; bk:=k; bS:=SetToSequence(T); bd:=e; end if;
    end if; end for; end for;
    return bk,bS,bd;
end function;
function Reps(F)
    A:={x:x in F|x ne F!0 and x ne F!1}; seen:=AssociativeArray(); R:=[];
    for T in Subsets(A,4) do for d in F do
        k,S,e:=CR(F,{F!0,F!1} join T,d);
        if not IsDefined(seen,k) then seen[k]:=true; Append(~R,<S,e>); end if;
    end for; end for; return R;
end function;
function RV(S,d)
    return #[I:I in Subsets({1..6},3)|&+[S[i]:i in SetToSequence(I)] eq d];
end function;
function NV(P)
    F:=BaseRing(Parent(P));
    for i in [1..4] do if P[i] ne 0 then return Vector(F,[P[j]/P[i]:j in [1..4]]); end if; end for;
    return P;
end function;
function PK(P) return Join([Sprint(x):x in Eltseq(NV(P))],","); end function;
procedure AS(~A,F,C)
    for u in CartesianPower(F,#C) do
        P:=Vector(F,[F!0,F!0,F!0,F!0]);
        for i in [1..#C] do P+:=u[i]*C[i]; end for;
        if not IsZero(P) then A[PK(P)]:=true; end if;
    end for;
end procedure;
procedure Scan(q)
    F<w>:=GF(q); R:=Reps(F); mism:=0; byR:=AssociativeArray();
    for X in R do
        S:=X[1]; d:=X[2]; th:=&+S-d; r:=RV(S,d);
        C:=[Vector(F,[F!1,a,a^2,a^3]):a in S] cat
          [Vector(F,[0,0,1,th]),Vector(F,[0,0,0,1])];
        A1:=AssociativeArray(); A2:=AssociativeArray();
        for i in [1..8] do AS(~A1,F,[C[i]]); end for;
        for I in Subsets({1..8},2) do AS(~A2,F,[C[i]:i in SetToSequence(I)]); end for;
        m2:=#SetToSequence(Keys(A2)) - #SetToSequence(Keys(A1));
        expected:=28*(q-1)-3*r;
        if m2 ne expected then mism+:=1; end if;
        k:=Sprint(r); if not IsDefined(byR,k) then byR[k]:=0; end if; byR[k]+:=1;
    end for;
    printf "M2_FORMULA q=%o classes=%o mismatches=%o rDistribution=%o\n",q,#R,mism,byR;
end procedure;
for q in [16,17,19,23] do Scan(q); end for;
quit;
