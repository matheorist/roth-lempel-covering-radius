// codex_rounds11_18_boundary_certificate_analysis.m
// Boundary analysis for the twisted-cubic certificate and the shifted family
// P(x,mu)=V(x)+mu F=(1,x,x^2,x^3+mu).

SetEchoInput(false);

function ElemKey(x)
    return Sprint(x);
end function;

function SeqKey(S)
    return Join([ElemKey(x) : x in S],",");
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
    return #[I : I in Subsets({1..#S},3) |
        &+[S[i] : i in SetToSequence(I)] eq delta];
end function;

function ForbiddenCubicParameters(S, delta)
    theta := &+[a : a in S]-delta;
    bad := SequenceToSet(S);
    for I in Subsets({1..#S},2) do
        J := SetToSequence(I);
        Include(~bad,theta-S[J[1]]-S[J[2]]);
    end for;
    return bad;
end function;

function ProductFactorsNonzero(P, S, delta)
    theta := &+[a : a in S]-delta;
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

function ShiftedCurveWitness(F, S, delta)
    for x in F do for mu in F do
        if ProductFactorsNonzero(Vector(F,[F!1,x,x^2,x^3+mu]),S,delta) then
            return true,x,mu;
        end if;
    end for; end for;
    return false,F!0,F!0;
end function;

function FullWitness(F, S, delta)
    for a in F do for b in F do for c in F do
        if ProductFactorsNonzero(Vector(F,[F!1,a,b,c]),S,delta) then
            return true,Vector(F,[F!1,a,b,c]);
        end if;
    end for; end for; end for;
    for a in F do for b in F do
        if ProductFactorsNonzero(Vector(F,[F!0,F!1,a,b]),S,delta) then
            return true,Vector(F,[F!0,F!1,a,b]);
        end if;
    end for; end for;
    for a in F do
        if ProductFactorsNonzero(Vector(F,[F!0,F!0,F!1,a]),S,delta) then
            return true,Vector(F,[F!0,F!0,F!1,a]);
        end if;
    end for;
    return false,Vector(F,[F!0,F!0,F!0,F!1]);
end function;

procedure BoundaryScan(q)
    F<w> := GF(q); reps := AffineRepresentatives(F);
    Bdist := AssociativeArray(); curveFails := []; shiftedFails := [];
    shiftedAdd := 0;
    for rep in reps do
        S := rep[1]; d := rep[2];
        r := TripleSumMultiplicity(S,d); bsz := #ForbiddenCubicParameters(S,d);
        key := "r=" cat Sprint(r) cat "|B=" cat Sprint(bsz);
        if not IsDefined(Bdist,key) then Bdist[key] := 0; end if;
        Bdist[key] +:= 1;
        if bsz eq q then
            Append(~curveFails,<S,d,r,bsz>);
            ok,x,mu := ShiftedCurveWitness(F,S,d);
            if ok then
                shiftedAdd +:= 1;
                printf "SHIFT_CERT q=%o S=%o delta=%o r=%o x=%o mu=%o\n",
                    q,S,d,r,x,mu;
            else
                ok2,P := FullWitness(F,S,d);
                Append(~shiftedFails,<S,d,r,P>);
                printf "FULL_ONLY q=%o S=%o delta=%o r=%o found=%o P=%o\n",
                    q,S,d,r,ok2,P;
            end if;
        end if;
    end for;
    printf "BOUNDARY_SUMMARY q=%o classes=%o curveFails=%o shiftedAdds=%o shiftedFails=%o\n",
        q,#reps,#curveFails,shiftedAdd,#shiftedFails;
    printf "B_DISTRIBUTION q=%o\n",q;
    for key in Sort(SetToSequence(Keys(Bdist))) do
        printf "%o classes=%o\n",key,Bdist[key];
    end for;
end procedure;

for q in [16,17,19,23] do BoundaryScan(q); end for;
quit;
