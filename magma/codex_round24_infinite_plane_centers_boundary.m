// codex_round24_infinite_plane_centers_boundary.m
// Tests the plane at infinity family Q(t,u)=(0,1,t,u) on curve failures.

SetEchoInput(false);
function EK(x) return Sprint(x); end function;
function CK(S,d) return Join(Sort([EK(x):x in S]),",") cat "|" cat EK(d); end function;
function CR(F,S,d)
    first:=true; k0:=""; S0:=S; d0:=d;
    for x in S do for y in S do if x ne y then
        a:=(y-x)^-1; b:=-a*x; T:={a*z+b:z in S}; e:=a*d+3*b; k:=CK(T,e);
        if first or k lt k0 then first:=false;k0:=k;S0:=SetToSequence(T);d0:=e;end if;
    end if; end for; end for; return k0,S0,d0;
end function;
function Reps(F)
    A:={x:x in F|x ne F!0 and x ne F!1}; seen:=AssociativeArray(); R:=[];
    for X in Subsets(A,4) do for d in F do
        k,S,e:=CR(F,{F!0,F!1} join X,d);
        if not IsDefined(seen,k) then seen[k]:=true;Append(~R,<S,e>);end if;
    end for;end for;return R;
end function;
function RVal(S,d) return #[I:I in Subsets({1..6},3)|&+[S[i]:i in SetToSequence(I)] eq d]; end function;
function BSet(S,d)
    th:=&+S-d; B:=SequenceToSet(S);
    for I in Subsets({1..6},2) do J:=SetToSequence(I);Include(~B,th-S[J[1]]-S[J[2]]);end for; return B;
end function;
function Hole(P,S,d)
    th:=&+S-d;
    for a in S do if P[2]-a*P[1] eq 0 then return false; end if;end for;
    for I in Subsets({1..6},2) do J:=SetToSequence(I);a:=S[J[1]];b:=S[J[2]];c:=th-a-b;
        if P[3]-(a+b)*P[2]+a*b*P[1] eq 0 then return false;end if;
        if P[4]-th*P[3]+(a*b+(a+b)*c)*P[2]-a*b*c*P[1] eq 0 then return false;end if;
    end for;
    for I in Subsets({1..6},3) do J:=SetToSequence(I);a:=S[J[1]];b:=S[J[2]];c:=S[J[3]];
        if P[4]-(a+b+c)*P[3]+(a*b+a*c+b*c)*P[2]-a*b*c*P[1] eq 0 then return false; end if;
    end for; return true;
end function;
procedure Scan(q)
    F<w>:=GF(q); failCurve:=0; Qadd:=0; Qfail:=0;
    for X in Reps(F) do S:=X[1];d:=X[2];
        if #BSet(S,d) eq q then
            failCurve+:=1; found:=false; tw:=F!0; uw:=F!0;
            for t in F do for u in F do
                if not found and Hole(Vector(F,[0,1,t,u]),S,d) then found:=true;tw:=t;uw:=u;end if;
            end for;end for;
            if found then Qadd+:=1; printf "Q_CERT q=%o S=%o delta=%o r=%o t=%o u=%o\n",q,S,d,RVal(S,d),tw,uw;
            else Qfail+:=1;end if;
        end if;
    end for;
    printf "Q_SUMMARY q=%o curveFails=%o Qcertified=%o Qfailed=%o\n",q,failCurve,Qadd,Qfail;
end procedure;
for q in [16,17,19,23] do Scan(q); end for;
quit;
