(* Dilithium Parameter Definitions *)
Require Import ZArith.

Module DilithiumParams.
  Definition N : Z := 256.
  Definition Q : Z := 8380417.
  Definition D : Z := 13.
  Definition ROOT_OF_UNITY : Z := 1753.

  (* Mode: mode3 *)
  Definition K : Z := 6.
  Definition L : Z := 5.
  Definition ETA : Z := 4.
  Definition TAU : Z := 49.
  Definition BETA : Z := 196.
  Definition GAMMA1 : Z := 524288.
  Definition GAMMA2 : Z := 261888.
  Definition OMEGA : Z := 55.

  (* Derived Sizes *)
  Definition CRYPTO_PUBLICKEYBYTES : Z := 1952.
  Definition CRYPTO_SECRETKEYBYTES : Z := 4032.
  Definition CRYPTO_BYTES : Z := 3309.
End DilithiumParams.
