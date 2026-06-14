// codex_phase9_lowfield_triplesum_hole_spectrum_codex.m
// Fast structural spectrum: number of three-subset representations of delta
// versus existence of distance-four syndrome cosets.

SetEchoInput(false);

function ElemKey(x)
    return Sprint(x);
end function;

function CaseKey(S, delta)
    return Join(Sort([ElemKey(x) : x in S]),",") cat "|" cat ElemKey(delta);
end function;

function CanonicalRepresentative(F, S, delta)
    first := true; bestKey := ""; bestS := S; bestD := delta;
    for x in S do for y in S do
        if x eq y then continue; end if;
        a := (y-x)^-1; b := -a*x;
        T := {a*z+b : z in S}; d := a*delta+3*b;
        key := CaseKey(T,d);
        if first or key lt bestKey then
            first := false; bestKey := key; bestS := SetToSequence(T); bestD := d;
        end if;
    end for; end for;
    return bestKey,bestS,bestD;
end function;

function AffineRepresentatives(F)
    base := {x : x in F | x ne F!0 and x ne F!1};
    seen := AssociativeArray(); reps := [];
    for tail in Subsets(base,4) do
        for delta in F do
            key,S,d := CanonicalRepresentative(F,{F!0,F!1} join tail,delta);
            if not IsDefined(seen,key) then
                seen[key] := true; Append(~reps,<S,d>);
            end if;
        end for;
    end for;
    return reps;
end function;

function TripleSumMultiplicity(S, delta)
    r := 0;
    for I in Subsets({1..#S},3) do
        if &+[S[i] : i in SetToSequence(I)] eq delta then r +:= 1; end if;
    end for;
    return r;
end function;

function IsUncovered(P, S, delta)
    theta := &+[a : a in S] - delta;
    for a in S do
        if P[2]-a*P[1] eq 0 then return false; end if;
    end for;
    for I in Subsets({1..#S},2) do
        J := SetToSequence(I); a := S[J[1]]; b := S[J[2]];
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 then return false; end if;
        c := theta-a-b;
        if P[4]-theta*P[3]+(a*b+(a+b)*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    for I in Subsets({1..#S},3) do
        J := SetToSequence(I); a := S[J[1]]; b := S[J[2]]; c := S[J[3]];
        if P[4]-(a+b+c)*P[3]+(a*b+a*c+b*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function HasHole(F, S, delta)
    for a in F do for b in F do for c in F do
        if IsUncovered(Vector(F,[F!1,a,b,c]),S,delta) then return true; end if;
    end for; end for; end for;
    for a in F do for b in F do
        if IsUncovered(Vector(F,[F!0,F!1,a,b]),S,delta) then return true; end if;
    end for; end for;
    for a in F do
        if IsUncovered(Vector(F,[F!0,F!0,F!1,a]),S,delta) then return true; end if;
    end for;
    return IsUncovered(Vector(F,[F!0,F!0,F!0,F!1]),S,delta);
end function;

procedure ScanQ(q)
    F<w> := GF(q); reps := AffineRepresentatives(F);
    spectrum := AssociativeArray();
    for rep in reps do
        S := rep[1]; delta := rep[2];
        r := TripleSumMultiplicity(S,delta);
        rho := HasHole(F,S,delta) select 4 else 3;
        key := "r=" cat Sprint(r) cat "|rho=" cat Sprint(rho);
        if not IsDefined(spectrum,key) then spectrum[key] := 0; end if;
        spectrum[key] +:= 1;
    end for;
    printf "SPECTRUM q=%o classes=%o\n",q,#reps;
    for key in Sort(SetToSequence(Keys(spectrum))) do
        printf "%o classes=%o\n",key,spectrum[key];
    end for;
end procedure;

for q in [7,8,9,11,13] do ScanQ(q); end for;
quit;
