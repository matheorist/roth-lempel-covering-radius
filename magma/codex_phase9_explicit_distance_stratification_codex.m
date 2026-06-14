// codex_phase9_explicit_distance_stratification_codex.m
// Full explicit syndrome-distance stratification and the MDS/NMDS
// triple-sum criterion for RL_{4,delta}(S), |S|=6.

SetEchoInput(false);

function ConstructRLCode(F, S, delta)
    rows := [];
    Append(~rows, [F | F!1 : a in S] cat [F!0,F!0]);
    Append(~rows, [F | a : a in S] cat [F!0,F!0]);
    Append(~rows, [F | a^2 : a in S] cat [F!0,F!1]);
    Append(~rows, [F | a^3 : a in S] cat [F!1,delta]);
    return LinearCode(Matrix(F,4,8,rows));
end function;

function ExplicitParityCheck(F, S, delta)
    sigma := &+[a : a in S];
    lambda := [(&*[S[i]-S[j] : j in [1..#S] | j ne i])^-1
        : i in [1..#S]];
    rows := [
        [F | lambda[i] : i in [1..#S]] cat [F!0,F!0],
        [F | lambda[i]*S[i] : i in [1..#S]] cat [F!0,F!0],
        [F | lambda[i]*S[i]^2 : i in [1..#S]] cat [-F!1,F!0],
        [F | lambda[i]*S[i]^3 : i in [1..#S]]
            cat [delta-sigma,-F!1]
    ];
    return Matrix(F,4,8,rows);
end function;

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
        alpha := (y-x)^-1; beta := -alpha*x;
        T := {alpha*z+beta : z in S};
        d := alpha*delta+3*beta;
        key := CaseKey(T,d);
        if first or key lt bestKey then
            first := false; bestKey := key; bestS := SetToSequence(T); bestD := d;
        end if;
    end for; end for;
    return bestKey, bestS, bestD;
end function;

function AffineRepresentatives(F)
    base := {x : x in F | x ne F!0 and x ne F!1};
    seen := AssociativeArray(); reps := [];
    for tail in Subsets(base,4) do
        S := {F!0,F!1} join tail;
        for delta in F do
            key,T,d := CanonicalRepresentative(F,S,delta);
            if not IsDefined(seen,key) then
                seen[key] := true; Append(~reps,<T,d>);
            end if;
        end for;
    end for;
    return reps;
end function;

function Col(M, c)
    F := BaseRing(M);
    return Vector(F,[M[r][c] : r in [1..Nrows(M)]]);
end function;

function ProjectivePoints3(F)
    pts := [Vector(F,[F!1,a,b,c]) : a,b,c in F];
    pts cat:= [Vector(F,[F!0,F!1,a,b]) : a,b in F];
    pts cat:= [Vector(F,[F!0,F!0,F!1,a]) : a in F];
    Append(~pts,Vector(F,[F!0,F!0,F!0,F!1]));
    return pts;
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
    return Join([Sprint(x) : x in Eltseq(NormalizeVector(P))],",");
end function;

procedure AddProjectiveSpan(~covered, F, cols)
    for coeff in CartesianPower(F,#cols) do
        P := Vector(F,[F!0 : i in [1..4]]);
        for i in [1..#cols] do P +:= coeff[i]*cols[i]; end for;
        if not IsZero(P) then covered[PointKey(P)] := true; end if;
    end for;
end procedure;

function GeometricShellMaps(H)
    F := BaseRing(H);
    cols := [Col(H,i) : i in [1..Ncols(H)]];
    maps := [AssociativeArray(), AssociativeArray(), AssociativeArray()];
    for t in [1..3] do
        for I in Subsets({1..#cols},t) do
            AddProjectiveSpan(~maps[t],F,[cols[i] : i in SetToSequence(I)]);
        end for;
    end for;
    return maps;
end function;

function IsNMDSByTripleSum(S, delta)
    for T in Subsets({1..#S},3) do
        I := SetToSequence(T);
        if &+[S[i] : i in I] eq delta then return true; end if;
    end for;
    return false;
end function;

function HasDistanceAtMostTwo(P, S, delta)
    theta := &+[a : a in S] - delta;
    if P[1] eq 0 and P[2] eq 0 then return true; end if; // <E,F>
    for a in S do
        if P[2]-a*P[1] eq 0 and P[3]-a^2*P[1] eq 0 then
            return true; // <V(a),F>
        end if;
        if P[2]-a*P[1] eq 0 and
            P[4]-theta*P[3]-(a^3-theta*a^2)*P[1] eq 0 then
            return true; // <V(a),E>
        end if;
    end for;
    for I in Subsets({1..#S},2) do
        inds := SetToSequence(I); a := S[inds[1]]; b := S[inds[2]];
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 and
            P[4]-(a+b)*P[3]+a*b*P[2] eq 0 then
            return true; // <V(a),V(b)>
        end if;
    end for;
    return false;
end function;

function IsUncovered(P, S, delta)
    theta := &+[a : a in S] - delta;
    for a in S do
        if P[2]-a*P[1] eq 0 then return false; end if;
    end for;
    for I in Subsets({1..#S},2) do
        inds := SetToSequence(I); a := S[inds[1]]; b := S[inds[2]];
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 then return false; end if;
        c := theta-a-b;
        if P[4]-theta*P[3]+(a*b+(a+b)*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    for I in Subsets({1..#S},3) do
        inds := SetToSequence(I);
        a := S[inds[1]]; b := S[inds[2]]; c := S[inds[3]];
        if P[4]-(a+b+c)*P[3]+(a*b+a*c+b*c)*P[2]-a*b*c*P[1] eq 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function HasDistanceAtMostOne(P, H)
    F := BaseRing(H);
    for i in [1..Ncols(H)] do
        if Rank(Matrix(F,2,4,Eltseq(P) cat Eltseq(Col(H,i)))) eq 1 then
            return true;
        end if;
    end for;
    return false;
end function;

function ExplicitDistance(P, H, S, delta)
    if HasDistanceAtMostOne(P,H) then return 1; end if;
    if HasDistanceAtMostTwo(P,S,delta) then return 2; end if;
    if IsUncovered(P,S,delta) then return 4; end if;
    return 3;
end function;

function GeometricDistance(P, maps)
    key := PointKey(P);
    for t in [1..3] do
        if IsDefined(maps[t],key) then return t; end if;
    end for;
    return 4;
end function;

procedure VerifyQ(q)
    F<w> := GF(q);
    reps := AffineRepresentatives(F); pts := ProjectivePoints3(F);
    typeMismatch := 0; dualCriterionMismatch := 0; distanceMismatch := 0;
    rho2 := 0; rho3 := 0; rho4 := 0;
    for rep in reps do
        S := rep[1]; delta := rep[2]; C := ConstructRLCode(F,S,delta);
        bySum := IsNMDSByTripleSum(S,delta);
        if bySum eq IsMDS(C) then typeMismatch +:= 1; end if;
        theta := &+[a : a in S]-delta;
        if IsNMDSByTripleSum(S,delta) ne IsNMDSByTripleSum(S,theta) then
            dualCriterionMismatch +:= 1;
        end if;
        H := ExplicitParityCheck(F,S,delta);
        maps := GeometricShellMaps(H);
        maxd := 0;
        for P in pts do
            de := ExplicitDistance(P,H,S,delta);
            dg := GeometricDistance(P,maps);
            if de ne dg then distanceMismatch +:= 1; end if;
            if de gt maxd then maxd := de; end if;
        end for;
        if maxd eq 2 then rho2 +:= 1;
        elif maxd eq 3 then rho3 +:= 1;
        elif maxd eq 4 then rho4 +:= 1;
        end if;
    end for;
    printf "SUMMARY q=%o classes=%o typeMismatch=%o dualCriterionMismatch=%o distanceMismatch=%o rho2=%o rho3=%o rho4=%o\n",
        q,#reps,typeMismatch,dualCriterionMismatch,distanceMismatch,rho2,rho3,rho4;
end procedure;

if not assigned runQs then runQs := [7,8,9,11,13,16]; end if;
for q in runQs do VerifyQ(q); end for;
quit;
