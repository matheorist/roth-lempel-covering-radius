// codex_round16_lowfield_explicit_shell_by_r_codex.m
// Low-field shell spectrum grouped by triple-sum multiplicity r_S(delta).

SetEchoInput(false);

function ElemKey(x)
    return Sprint(x);
end function;

function CaseKey(S, d)
    return Join(Sort([ElemKey(x) : x in S]),",") cat "|" cat ElemKey(d);
end function;

function CanonicalRepresentative(F, S, d)
    first := true; bestKey := ""; bestS := S; bestD := d;
    for x in S do for y in S do
        if x eq y then continue; end if;
        a := (y-x)^-1; b := -a*x;
        T := {a*z+b : z in S}; e := a*d+3*b;
        key := CaseKey(T,e);
        if first or key lt bestKey then
            first := false; bestKey := key; bestS := SetToSequence(T); bestD := e;
        end if;
    end for; end for;
    return bestKey,bestS,bestD;
end function;

function AffineRepresentatives(F)
    base := {x : x in F | x ne F!0 and x ne F!1};
    seen := AssociativeArray(); reps := [];
    for tail in Subsets(base,4) do for d in F do
        key,S,e := CanonicalRepresentative(F,{F!0,F!1} join tail,d);
        if not IsDefined(seen,key) then seen[key] := true; Append(~reps,<S,e>); end if;
    end for; end for;
    return reps;
end function;

function RValue(S,d)
    return #[I : I in Subsets({1..6},3) |
        &+[S[i] : i in SetToSequence(I)] eq d];
end function;

function NormalizedColumns(F,S,d)
    theta := &+S-d;
    return [Vector(F,[F!1,a,a^2,a^3]) : a in S] cat
        [Vector(F,[F!0,F!0,F!1,theta]),Vector(F,[F!0,F!0,F!0,F!1])];
end function;

function Normalize(P)
    F := BaseRing(Parent(P));
    for i in [1..4] do
        if P[i] ne 0 then return Vector(F,[P[j]/P[i] : j in [1..4]]); end if;
    end for;
    return P;
end function;

function Key(P)
    return Join([Sprint(x) : x in Eltseq(Normalize(P))],",");
end function;

procedure AddSpan(~A,F,cols)
    for coeff in CartesianPower(F,#cols) do
        P := Vector(F,[F!0,F!0,F!0,F!0]);
        for i in [1..#cols] do P +:= coeff[i]*cols[i]; end for;
        if not IsZero(P) then A[Key(P)] := true; end if;
    end for;
end procedure;

function Points(F)
    P := [Vector(F,[F!1,a,b,c]) : a,b,c in F];
    P cat:= [Vector(F,[F!0,F!1,a,b]) : a,b in F];
    P cat:= [Vector(F,[F!0,F!0,F!1,a]) : a in F];
    Append(~P,Vector(F,[F!0,F!0,F!0,F!1]));
    return P;
end function;

function IsHole(P,S,d)
    theta := &+S-d;
    for a in S do if P[2]-a*P[1] eq 0 then return false; end if; end for;
    for I in Subsets({1..6},2) do
        J := SetToSequence(I); a:=S[J[1]]; b:=S[J[2]]; c:=theta-a-b;
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 then return false; end if;
        if P[4]-theta*P[3]+(a*b+(a+b)*c)*P[2]-a*b*c*P[1] eq 0 then return false; end if;
    end for;
    for I in Subsets({1..6},3) do
        J := SetToSequence(I); a:=S[J[1]]; b:=S[J[2]]; c:=S[J[3]];
        if P[4]-(a+b+c)*P[3]+(a*b+a*c+b*c)*P[2]-a*b*c*P[1] eq 0 then return false; end if;
    end for;
    return true;
end function;

procedure Scan(q)
    F<w> := GF(q); reps:=AffineRepresentatives(F); pts:=Points(F);
    group := AssociativeArray();
    for rep in reps do
        S:=rep[1]; d:=rep[2]; r:=RValue(S,d); cols:=NormalizedColumns(F,S,d);
        one:=AssociativeArray(); two:=AssociativeArray();
        for i in [1..8] do AddSpan(~one,F,[cols[i]]); end for;
        for I in Subsets({1..8},2) do AddSpan(~two,F,[cols[i] : i in SetToSequence(I)]); end for;
        shells:=[0,0,0,0];
        for P in pts do
            if IsDefined(one,Key(P)) then shells[1]+:=1;
            elif IsDefined(two,Key(P)) then shells[2]+:=1;
            elif IsHole(P,S,d) then shells[4]+:=1;
            else shells[3]+:=1;
            end if;
        end for;
        rho := shells[4] gt 0 select 4 else 3;
        key := "r=" cat Sprint(r) cat "|rho=" cat Sprint(rho) cat "|shells=" cat Sprint(shells);
        if not IsDefined(group,key) then group[key]:=0; end if;
        group[key]+:=1;
    end for;
    printf "SHELL_SPECTRUM q=%o classes=%o\n",q,#reps;
    for key in Sort(SetToSequence(Keys(group))) do printf "%o classes=%o\n",key,group[key]; end for;
end procedure;

if not assigned runQs then runQs := [7,8,9,11,13]; end if;
for q in runQs do Scan(q); end for;
quit;
