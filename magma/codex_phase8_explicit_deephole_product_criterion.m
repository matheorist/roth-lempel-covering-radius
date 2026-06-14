// codex_phase8_explicit_deephole_product_criterion.m
// Verification of an explicit parity-check matrix and the 56-factor
// uncovered-syndrome criterion for RL_{4,delta}(S), |S|=6.

SetEchoInput(false);

function ConstructRLGenerator(F, S, delta)
    rows := [];
    Append(~rows, [F | F!1 : a in S] cat [F!0,F!0]);
    Append(~rows, [F | a : a in S] cat [F!0,F!0]);
    Append(~rows, [F | a^2 : a in S] cat [F!0,F!1]);
    Append(~rows, [F | a^3 : a in S] cat [F!1,delta]);
    return Matrix(F,4,8,rows);
end function;

function ExplicitParityCheck(F, S, delta)
    sigma := &+[a : a in S];
    lambda := [];
    for i in [1..#S] do
        Append(~lambda, (&*[S[i]-S[j] : j in [1..#S] | j ne i])^-1);
    end for;
    rows := [];
    Append(~rows, [F | lambda[i] : i in [1..#S]] cat [F!0,F!0]);
    Append(~rows, [F | lambda[i]*S[i] : i in [1..#S]] cat [F!0,F!0]);
    Append(~rows, [F | lambda[i]*S[i]^2 : i in [1..#S]] cat [-F!1,F!0]);
    Append(~rows, [F | lambda[i]*S[i]^3 : i in [1..#S]]
        cat [delta-sigma,-F!1]);
    return Matrix(F,4,8,rows);
end function;

function ElemKey(x)
    return Sprint(x);
end function;

function SetKey(S)
    return Join(Sort([ElemKey(x) : x in S]),",");
end function;

function CaseKey(S, delta)
    return SetKey(S) cat "|" cat ElemKey(delta);
end function;

function CanonicalRepresentative(F, S, delta)
    first := true; bestKey := ""; bestS := S; bestD := delta;
    for x in S do
        for y in S do
            if x eq y then continue; end if;
            a := (y-x)^-1; b := -a*x;
            T := {a*z+b : z in S};
            d := a*delta+3*b;
            key := CaseKey(T,d);
            if first or key lt bestKey then
                first := false; bestKey := key;
                bestS := SetToSequence(T); bestD := d;
            end if;
        end for;
    end for;
    return bestKey, bestS, bestD;
end function;

function AffineRepresentatives(F)
    zero := F!0; one := F!1;
    base := {x : x in F | x ne zero and x ne one};
    seen := AssociativeArray(); reps := [];
    for tail in Subsets(base,4) do
        S := {zero,one} join tail;
        for delta in F do
            key,T,d := CanonicalRepresentative(F,S,delta);
            if not IsDefined(seen,key) then
                seen[key] := true;
                Append(~reps,<T,d>);
            end if;
        end for;
    end for;
    return reps;
end function;

function Col(M, c)
    F := BaseRing(M);
    return Vector(F, [M[r][c] : r in [1..Nrows(M)]]);
end function;

function TriplePlaneEquations(H)
    F := BaseRing(H);
    cols := [Col(H,i) : i in [1..Ncols(H)]];
    eqns := [ VectorSpace(F,4) | ];
    for ss in Subsets({1..Ncols(H)},3) do
        tri := Sort(SetToSequence(ss));
        A := Matrix(F,3,4,Eltseq(cols[tri[1]]) cat
            Eltseq(cols[tri[2]]) cat Eltseq(cols[tri[3]]));
        if Rank(A) lt 3 then return false, eqns; end if;
        Append(~eqns, Basis(Nullspace(Transpose(A)))[1]);
    end for;
    return true, eqns;
end function;

function IsGeometricallyUncovered(P, eqns)
    for e in eqns do
        if &+[e[i]*P[i] : i in [1..4]] eq 0 then return false; end if;
    end for;
    return true;
end function;

function IsExplicitlyUncovered(P, S, delta)
    theta := &+[a : a in S] - delta;
    for a in S do
        if P[2]-a*P[1] eq 0 then return false; end if;
    end for;
    for ij in Subsets({1..#S},2) do
        indices := Sort(SetToSequence(ij));
        a := S[indices[1]]; b := S[indices[2]];
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 then return false; end if;
        c := theta-a-b;
        if P[4]-theta*P[3]+(a*b+(a+b)*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    for ijk in Subsets({1..#S},3) do
        indices := Sort(SetToSequence(ijk));
        a := S[indices[1]]; b := S[indices[2]]; c := S[indices[3]];
        if P[4]-(a+b+c)*P[3]+(a*b+a*c+b*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function ProjectivePoints3(F)
    pts := [ Vector(F,[F!1,a,b,c]) : a,b,c in F ];
    pts cat:= [ Vector(F,[F!0,F!1,a,b]) : a,b in F ];
    pts cat:= [ Vector(F,[F!0,F!0,F!1,a]) : a in F ];
    Append(~pts, Vector(F,[F!0,F!0,F!0,F!1]));
    return pts;
end function;

procedure VerifyQ(q)
    F<w> := GF(q);
    reps := AffineRepresentatives(F);
    pts := ProjectivePoints3(F);
    badOrthogonal := 0; badRank := 0; badTriples := 0; mismatches := 0;
    totalUncovered := 0;
    printf "=== Explicit deep-hole product verification q=%o generator=%o ===\n", q, w;
    printf "AFFINE_CLASSES=%o PG_POINTS=%o\n", #reps, #pts;
    for rep in reps do
        S := rep[1]; delta := rep[2];
        G := ConstructRLGenerator(F,S,delta);
        H := ExplicitParityCheck(F,S,delta);
        if G*Transpose(H) ne ZeroMatrix(F,4,4) then badOrthogonal +:= 1; end if;
        if Rank(H) ne 4 then badRank +:= 1; end if;
        valid, eqns := TriplePlaneEquations(H);
        if not valid then
            badTriples +:= 1;
        else
            for P in pts do
                geom := IsGeometricallyUncovered(P,eqns);
                expl := IsExplicitlyUncovered(P,S,delta);
                if geom then totalUncovered +:= 1; end if;
                if geom ne expl then mismatches +:= 1; end if;
            end for;
        end if;
    end for;
    printf "SUMMARY q=%o badOrthogonal=%o badRank=%o badTriples=%o mismatches=%o totalUncoveredOverClasses=%o\n",
        q,badOrthogonal,badRank,badTriples,mismatches,totalUncovered;
end procedure;

if not assigned runQs then runQs := [7,8,9,11,13,16]; end if;
for q in runQs do VerifyQ(q); end for;
quit;
