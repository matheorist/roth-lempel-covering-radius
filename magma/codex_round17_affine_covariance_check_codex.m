// codex_round17_affine_covariance_check_codex.m
// Checks covariance of theta and the twisted-cubic forbidden set.

SetEchoInput(false);

function Forbidden(S,d)
    th:=&+S-d; A:=SequenceToSet(S);
    for I in Subsets({1..#S},2) do
        J:=SetToSequence(I); Include(~A,th-S[J[1]]-S[J[2]]);
    end for;
    return A;
end function;

procedure Check(q)
    F<w>:=GF(q); failTheta:=0; failForbidden:=0; trials:=0;
    seq:=SetToSequence(Set(F));
    for trial in [1..50] do
        S:=[seq[i] : i in [1..6]];
        // Change the starting six-set deterministically across trials.
        S:=[seq[((i+trial-2) mod q)+1] : i in [1..6]];
        if #SequenceToSet(S) lt 6 then continue; end if;
        d:=seq[((trial*3) mod q)+1];
        for a in F do if a ne 0 then
            b:=seq[((trial+Index(seq,a)) mod q)+1];
            T:=[a*x+b:x in S]; e:=a*d+3*b;
            th:=&+S-d; th2:=&+T-e;
            if th2 ne a*th+3*b then failTheta+:=1; end if;
            image:={a*x+b:x in Forbidden(S,d)};
            if image ne Forbidden(T,e) then failForbidden+:=1; end if;
            trials+:=1;
        end if; end for;
    end for;
    printf "AFFINE_COVARIANCE q=%o trials=%o thetaFailures=%o forbiddenFailures=%o\n",
        q,trials,failTheta,failForbidden;
end procedure;

for q in [7,8,9,11,13,16,17,19,23] do Check(q); end for;
quit;
