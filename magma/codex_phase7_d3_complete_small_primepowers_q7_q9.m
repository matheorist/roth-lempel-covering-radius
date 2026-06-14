// codex_phase7_d3_complete_small_primepowers_q7_q9.m
// Exact affine classification for the missing admissible fields q=7,8,9.
// For these small fields we calculate CoveringRadius for every class, rather
// than inferring only rho=4 from an uncovered syndrome certificate.

SetEchoInput(false);

function ConstructRLCode(F, S, k, delta)
    rows := [];
    for r in [0..k-1] do
        row := [F | S[i]^r : i in [1..#S]];
        if r eq k-1 then
            Append(~row,F!1); Append(~row,delta);
        elif r eq k-2 then
            Append(~row,F!0); Append(~row,F!1);
        else
            Append(~row,F!0); Append(~row,F!0);
        end if;
        Append(~rows,row);
    end for;
    return LinearCode(Matrix(F,k,#S+2,rows));
end function;

function ElemKey(x)
    return Sprint(x);
end function;

function SeqKey(seq)
    return Join([ElemKey(x) : x in seq],",");
end function;

function SetKey(S)
    return Join(Sort([ElemKey(x) : x in S]),",");
end function;

function CaseKey(S, delta)
    return SetKey(S) cat "|" cat ElemKey(delta);
end function;

function CanonicalRepresentative(F, S, delta, k)
    first := true; bestKey := ""; bestS := S; bestD := delta;
    for x in S do
        for y in S do
            if x eq y then continue; end if;
            a := (y-x)^-1; b := -a*x;
            T := {a*z+b : z in S};
            d := a*delta+(k-1)*b;
            key := CaseKey(T,d);
            if first or key lt bestKey then
                first := false; bestKey := key; bestS := SetToSequence(T); bestD := d;
            end if;
        end for;
    end for;
    return bestKey, bestS, bestD;
end function;

function AffineStabilizerSize(F, S, delta, k)
    count := 0;
    for a in F do
        if a eq 0 then continue; end if;
        for b in F do
            if {a*x+b : x in S} eq S and a*delta+(k-1)*b eq delta then
                count +:= 1;
            end if;
        end for;
    end for;
    return count;
end function;

function CodeType(C)
    N := Length(C); K := Dimension(C);
    if IsMDS(C) then return "MDS"; end if;
    if MinimumDistance(C) eq N-K and MinimumDistance(Dual(C)) eq K then
        return "NMDS";
    end if;
    return "OTHER";
end function;

function Col(M, c)
    F := BaseRing(M);
    return Vector(F, [M[r][c] : r in [1..Nrows(M)]]);
end function;

function NormalizeVector(P)
    F := BaseRing(Parent(P));
    for i in [1..#Eltseq(P)] do
        if P[i] ne 0 then
            return Vector(F,[P[j]/P[i] : j in [1..#Eltseq(P)]]);
        end if;
    end for;
    return P;
end function;

function PointKey(P)
    return SeqKey(Eltseq(NormalizeVector(P)));
end function;

procedure AddProjectiveSpan(~covered, F, cols)
    t := #cols;
    for coeff in CartesianPower(F,t) do
        P := Vector(F,[F!0 : i in [1..4]]);
        for i in [1..t] do P +:= coeff[i]*cols[i]; end for;
        if not IsZero(P) then covered[PointKey(P)] := true; end if;
    end for;
end procedure;

function DistanceShellCounts(C)
    F := BaseRing(C); H := GeneratorMatrix(Dual(C));
    cols := [Col(H,i) : i in [1..Ncols(H)]];
    covered1 := AssociativeArray(); covered2 := AssociativeArray(); covered3 := AssociativeArray();
    for ss in Subsets({1..Ncols(H)},1) do
        AddProjectiveSpan(~covered1,F,[cols[i] : i in SetToSequence(ss)]);
    end for;
    for ss in Subsets({1..Ncols(H)},2) do
        AddProjectiveSpan(~covered2,F,[cols[i] : i in SetToSequence(ss)]);
    end for;
    for ss in Subsets({1..Ncols(H)},3) do
        AddProjectiveSpan(~covered3,F,[cols[i] : i in SetToSequence(ss)]);
    end for;
    shells := [0,0,0,0];
    pts := [ Vector(F,[F!1,a,b,c]) : a,b,c in F ];
    pts cat:= [ Vector(F,[F!0,F!1,a,b]) : a,b in F ];
    pts cat:= [ Vector(F,[F!0,F!0,F!1,a]) : a in F ];
    Append(~pts,Vector(F,[F!0,F!0,F!0,F!1]));
    for P in pts do
        key := PointKey(P);
        if IsDefined(covered1,key) then shells[1] +:= 1;
        elif IsDefined(covered2,key) then shells[2] +:= 1;
        elif IsDefined(covered3,key) then shells[3] +:= 1;
        else shells[4] +:= 1;
        end if;
    end for;
    return shells;
end function;

function TriplePlaneEquations(C)
    F := BaseRing(C); H := GeneratorMatrix(Dual(C));
    cols := [Col(H,i) : i in [1..Ncols(H)]];
    eqns := [ VectorSpace(F,4) | ];
    for ss in Subsets({1..Ncols(H)},3) do
        tri := Sort(SetToSequence(ss));
        A := Matrix(F,3,4,Eltseq(cols[tri[1]]) cat Eltseq(cols[tri[2]]) cat Eltseq(cols[tri[3]]));
        if Rank(A) lt 3 then return false, eqns; end if;
        Append(~eqns, Basis(Nullspace(Transpose(A)))[1]);
    end for;
    return true, eqns;
end function;

function IsHole(P, eqns)
    for e in eqns do
        if &+[e[i]*P[i] : i in [1..4]] eq 0 then return false; end if;
    end for;
    return true;
end function;

function HoleCount(eqns, F)
    holes := 0;
    for a in F do for b in F do for c in F do
        if IsHole(Vector(F,[F!1,a,b,c]),eqns) then holes +:= 1; end if;
    end for; end for; end for;
    for a in F do for b in F do
        if IsHole(Vector(F,[F!0,F!1,a,b]),eqns) then holes +:= 1; end if;
    end for; end for;
    for a in F do
        if IsHole(Vector(F,[F!0,F!0,F!1,a]),eqns) then holes +:= 1; end if;
    end for;
    if IsHole(Vector(F,[F!0,F!0,F!0,F!1]),eqns) then holes +:= 1; end if;
    return holes;
end function;

procedure ScanQ(q)
    F<w> := GF(q); k := 4; n := 6;
    zero := F!0; one := F!1;
    base := {x : x in F | x ne zero and x ne one};
    seen := AssociativeArray(); reps := [];
    for tail in Subsets(base,4) do
        S := {zero,one} join tail;
        for delta in F do
            key,T,d := CanonicalRepresentative(F,S,delta,k);
            if not IsDefined(seen,key) then
                seen[key] := true;
                stab := AffineStabilizerSize(F,SequenceToSet(T),d,k);
                Append(~reps,<T,d,(q*(q-1)) div stab>);
            end if;
        end for;
    end for;

    rho2c := 0; rho2r := 0; rho3c := 0; rho3r := 0; rho4c := 0; rho4r := 0;
    mdsc := 0; nmdsc := 0; otherc := 0; invalidTriples := 0;
    spectrumClasses := AssociativeArray(); spectrumRaw := AssociativeArray();
    printf "=== Codex exact small-prime-power scan q=%o,k=4,n=6 generator=%o ===\n",q,w;
    printf "TOTAL_RAW=%o TOTAL_AFFINE_CLASSES=%o PG_POINTS=%o\n",
        Binomial(q,n)*q,#reps,q^3+q^2+q+1;
    for idx in [1..#reps] do
        S := reps[idx][1]; dlt := reps[idx][2]; orb := reps[idx][3];
        C := ConstructRLCode(F,S,k,dlt); typ := CodeType(C);
        if typ eq "MDS" then mdsc +:= 1;
        elif typ eq "NMDS" then nmdsc +:= 1;
        else otherc +:= 1;
        end if;
        valid,eqns := TriplePlaneEquations(C);
        if not valid then invalidTriples +:= 1; end if;
        holes := valid select HoleCount(eqns,F) else -1;
        rho := CoveringRadius(C);
        shells := DistanceShellCounts(C);
        deepShell := shells[rho];
        deepCosets := (q-1)*deepShell;
        deepWords := q^4*deepCosets;
        spkey := typ cat "|rho=" cat Sprint(rho) cat "|holes=" cat Sprint(holes);
        if not IsDefined(spectrumClasses,spkey) then
            spectrumClasses[spkey] := 0; spectrumRaw[spkey] := 0;
        end if;
        spectrumClasses[spkey] +:= 1; spectrumRaw[spkey] +:= orb;
        if rho eq 2 then rho2c +:= 1; rho2r +:= orb;
        elif rho eq 3 then rho3c +:= 1; rho3r +:= orb;
        elif rho eq 4 then rho4c +:= 1; rho4r +:= orb;
        end if;
        printf "CLASS idx=%o orbit=%o type=%o d=%o ddual=%o rho=%o holes=%o shells=%o deepCosets=%o deepWords=%o S=%o delta=%o\n",
            idx,orb,typ,MinimumDistance(C),MinimumDistance(Dual(C)),rho,holes,
            shells,deepCosets,deepWords,S,dlt;
    end for;
    printf "SUMMARY q=%o rho2Classes=%o rho2Raw=%o rho3Classes=%o rho3Raw=%o rho4Classes=%o rho4Raw=%o mdsClasses=%o nmdsClasses=%o otherClasses=%o invalidTriples=%o\n",
        q,rho2c,rho2r,rho3c,rho3r,rho4c,rho4r,mdsc,nmdsc,otherc,invalidTriples;
    printf "DEFECT_SPECTRUM q=%o\n",q;
    for key in Sort(SetToSequence(Keys(spectrumClasses))) do
        printf "%o classes=%o raw=%o\n",key,spectrumClasses[key],spectrumRaw[key];
    end for;
    printf "=== End q=%o exact small-prime-power scan ===\n",q;
end procedure;

if not assigned runQs then runQs := [7,8,9]; end if;
for q in runQs do ScanQ(q); end for;
quit;
