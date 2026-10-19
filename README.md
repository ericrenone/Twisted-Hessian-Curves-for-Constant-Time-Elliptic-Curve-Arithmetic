# Twisted Hessian Curves for Constant-Time Elliptic-Curve Arithmetic

A reference document for implementing the group law of a twisted Hessian elliptic curve in a hardware security module (HSM) or other constant-time environment. It gives the curve, the formulas, their operation counts, the conditions under which the formulas are complete, what the unified addition law does and does not buy against side channels, and a short test script that checks every formula in this document over a small prime field.

Everything here is either standard published material (cited) or a direct computation that the script in the appendix verifies. Nothing is claimed beyond that.

---

## Contents

1. [Scope](#1-scope)
2. [Notation](#2-notation)
3. [The curve](#3-the-curve)
4. [The group law](#4-the-group-law)
5. [Operation counts and comparison with other models](#5-operation-counts-and-comparison-with-other-models)
6. [Scalar multiplication](#6-scalar-multiplication)
7. [Side-channel properties](#7-side-channel-properties)
8. [Parameter requirements](#8-parameter-requirements)
9. [Implementation notes for an HSM](#9-implementation-notes-for-an-hsm)
10. [Known pitfalls](#10-known-pitfalls)
11. [Verification](#11-verification)
12. [Limitations and open items](#12-limitations-and-open-items)
13. [References](#13-references)
- [Appendix: verification script](#appendix-verification-script)

---

## 1. Scope

**What this covers**

- The projective twisted Hessian model `aX³ + Y³ + Z³ = dXYZ` over a prime field `F_p` with `p > 3`.
- The identity, negation, addition and doubling formulas, with operation counts.
- The condition under which a single addition formula is *complete* (works for every pair of input points, including doubling and the identity) so that a scalar multiplication can run the same instruction stream on every step.
- What that uniform instruction stream does and does not protect against.
- Parameter selection constraints specific to this model.

**What this does not cover**

- Field arithmetic itself (Montgomery multiplication, reduction, constant-time inversion). The document counts field operations and assumes they are constant time.
- Curve generation for production. The counting algorithm and security checks are named, not implemented.
- Any claim that this model is faster than twisted Edwards curves. It is not; see Section 5.
- Isogeny-based or post-quantum cryptography. Twisted Hessian curves are classical elliptic curves and the discrete logarithm on them is broken by a quantum computer like any other.

---

## 2. Notation

| Symbol | Meaning |
|---|---|
| `F_p` | prime field, `p > 3` |
| `[X : Y : Z]` | projective point, `(X, Y, Z) ≠ (0, 0, 0)`, equal up to a common non-zero scalar |
| `M` | one field multiplication |
| `S` | one field squaring |
| `Ma` | one multiplication by the curve constant `a` (free if `a = 1`, cheap if `a` is small) |
| `A` | one field addition or subtraction (not counted) |
| `O` | the group identity |
| `#E` | number of `F_p`-rational points on the curve |

Operation counts follow the convention of the Explicit-Formulas Database [EFD]: multiplications and squarings are counted, additions are not, and multiplications by small constants are listed separately.

---

## 3. The curve

### 3.1 Definition

For `a, d ∈ F_p`, the **twisted Hessian curve** is

```
TH(a,d):   aX³ + Y³ + Z³ = d·XYZ          (projective)
           ax³ + y³ + 1  = d·xy            (affine chart Z = 1, x = X/Z, y = Y/Z)
```

The curve is smooth (a genuine elliptic curve) exactly when

```
Δ(a,d) = a · (27a − d³) ≠ 0
```

that is, when `a ≠ 0` and `d³ ≠ 27a`. The sign convention for `Δ` varies between sources; only non-vanishing matters.

With `a = 1` this is the classical Hessian form `X³ + Y³ + Z³ = dXYZ` studied by Hesse (1844) and used in cryptography by Smart [Sma01] and Joye–Quisquater [JQ01]. The parameter `a` is the "twist": it lets the model cover curves whose extra flex points are not rational over `F_p`. Bernstein, Chuengsatiansup, Kohel and Lange [BCKL15] show that every elliptic curve over `F_p` with a rational point of order 3 is isomorphic to some `TH(a,d)`.

### 3.2 j-invariant

```
j(TH(a,d)) = d³ · (d³ + 216a)³ / ( a · (d³ − 27a)³ )
```

Sanity checks (both verified in the appendix script by a supersingularity test over `p ≡ 2 (mod 3)`):

- `d = 0` gives `j = 0` (the Fermat cubic `aX³ + Y³ + Z³ = 0`).
- `d³ = −216a` gives `j = 0`. With `a = 1`, that is `d = −6`.
- `d³ = 27a` is the singular locus (pole of `j`).

Note the **plus** sign in `(d³ + 216a)`. A version with a minus sign appears in some informal write-ups and is wrong.

### 3.3 Identity

```
O = [0 : −1 : 1]        (affine: (0, −1))
```

Check: `a·0 + (−1)³ + 1 = 0 = d·0·(−1)·1`.

`O` is an inflection point (flex), which is what makes the chord-and-tangent construction below produce a group with `O` as neutral element. Two independent checks:

*Hessian determinant.* For `F = aX³ + Y³ + Z³ − dXYZ` the matrix of second partials at `[0 : −1 : 1]` is

```
      F_XX  F_XY  F_XZ        0    −d    d
      F_YX  F_YY  F_YZ   =   −d    −6    0
      F_ZX  F_ZY  F_ZZ        d     0    6
```

(using `F_XX = 6aX`, `F_YY = 6Y`, `F_ZZ = 6Z`, `F_XY = −dZ`, `F_XZ = −dY`, `F_YZ = −dX`), whose determinant is `0·(−36) − (−d)(−6d − 0) + d(0 + 6d) = −6d² + 6d² = 0`. A point of a smooth cubic is a flex exactly when the Hessian determinant vanishes there.

*Tangent line.* The gradient at `O` is `(−d·(−1)·1, 3·1, 3·1) = (d, 3, 3)`, so the tangent line is `dX + 3Y + 3Z = 0`. Substituting a parametrisation of that line into `F` gives `(27a − d³)·t³` up to a non-zero constant, a triple root. Since `27a ≠ d³`, the tangent meets the curve only at `O`, with multiplicity 3.

### 3.4 Negation

```
−[X : Y : Z] = [X : Z : Y]            (affine: −(x, y) = (x/y, 1/y))
```

The curve equation is symmetric in `Y` and `Z`, so the swap preserves it. `O` is fixed by the swap, as it must be.

### 3.5 Flex points and 3-torsion

A smooth plane cubic over the algebraic closure has exactly nine flexes, and with a flex as identity they are exactly the points of order dividing 3. The rational ones on `TH(a,d)` are easy to locate:

- On the line `X = 0`: `Y³ + Z³ = 0`, so `Y/Z` is a cube root of `−1`. If `p ≡ 1 (mod 3)` there are three such roots and three rational points `[0 : −ω^i : 1]`, namely `O` and two points of order 3. If `p ≡ 2 (mod 3)` only `O` is rational.
- On the line `Z = 0`: `aX³ + Y³ = 0`, so `(Y/X)³ = −a`. These points are rational exactly when `−a` is a cube in `F_p`. When `a` is **not** a cube and `p ≡ 1 (mod 3)`, there are no rational points at infinity, which is the condition that makes the addition formula complete (Section 4.4).

Whatever the parameters, the group order satisfies

```
#TH(a,d)(F_p) ≡ 0  (mod 3)
```

This is a theorem in [BCKL15] and is reproduced numerically in the appendix for 24 parameter sets. Its practical consequence is stated in Section 8: the cofactor is at least 3, never 1.

---

## 4. The group law

Geometrically, three points of the curve sum to `O` exactly when they are collinear (with the usual multiplicity conventions), because `O` is a flex. So `P₁ + P₂` is the negation of the third intersection of the line `P₁P₂` with the curve. The formulas below are the explicit rational maps that compute this. All are taken from [BCKL15]; the operation counts are re-derived here and the formulas are tested in the appendix.

### 4.1 Standard addition

```
X₃ = X₁² Y₂ Z₂  −  X₂² Y₁ Z₁
Y₃ = Z₁² X₂ Y₂  −  Z₂² X₁ Y₁
Z₃ = Y₁² X₂ Z₂  −  Y₂² X₁ Z₁
```

**Cost: 12M.** Each of the six monomials is a product of two cross terms, e.g. `X₁²Y₂Z₂ = (X₁Y₂)(X₁Z₂)`. The six distinct cross terms `X₁Y₂, X₁Z₂, X₂Y₁, X₂Z₁, Y₂Z₁, Y₁Z₂` cost 6M, and combining them into the six monomials costs another 6M. No squarings, and the constant `a` does not appear.

**Exceptional cases.** This formula returns `(0, 0, 0)` when `P₁ = P₂` (so it cannot double), when either input has `X = 0` (so it cannot add `O`), and more generally when `P₁ − P₂` is a point of order 3. It is a fast *dedicated* addition, not a unified one.

### 4.2 Rotated addition (the complete formula)

```
X₃ = Z₂² X₁ Z₁  −  Y₁² X₂ Y₂
Y₃ = Y₂² Y₁ Z₁  −  a X₁² X₂ Z₂
Z₃ = a X₂² X₁ Y₁  −  Z₁² Y₂ Z₂
```

**Cost: 12M + 1Ma.** Factor each monomial into two cross terms:

```
Z₂²X₁Z₁  = (X₁Z₂)(Z₁Z₂)         Y₁²X₂Y₂  = (X₂Y₁)(Y₁Y₂)
Y₂²Y₁Z₁  = (Y₁Y₂)(Y₂Z₁)         aX₁²X₂Z₂ = (a·X₁X₂)(X₁Z₂)
aX₂²X₁Y₁ = (a·X₁X₂)(X₂Y₁)       Z₁²Y₂Z₂  = (Y₂Z₁)(Z₁Z₂)
```

Six distinct cross terms (`X₁X₂, Y₁Y₂, Z₁Z₂, X₁Z₂, X₂Y₁, Y₂Z₁`) cost 6M; `a·X₁X₂` costs 1Ma and is reused; the six final products cost 6M. With `a = 1` the formula is a flat 12M.

**Why it matters.** [BCKL15, Theorem 4.3] prove that this formula is **complete** whenever `a` is not a cube in `F_p`: for every pair of rational points, including `P₂ = P₁` (doubling), `P₂ = O`, and `P₂ = −P₁`, the output is a valid representative of `P₁ + P₂`, never `(0, 0, 0)`. The appendix confirms this on 40,000 pairs for a non-cube `a` and finds zero exceptions.

A complete addition law means one formula, one instruction stream, and no case analysis in the scalar-multiplication loop. That is the entire reason to prefer this model in a constant-time setting.

### 4.3 Dedicated doubling

Setting `P₂ = P₁` in the rotated formula and factoring:

```
X₃ = X₁ · (Z₁³ − Y₁³)
Y₃ = Z₁ · (Y₁³ − a X₁³)
Z₃ = Y₁ · (a X₁³ − Z₁³)
```

**Cost: 6M + 3S + 1Ma.** Three cubes (`1S + 1M` each) cost 3M + 3S, `a·X₁³` costs 1Ma, and the three outer products cost 3M.

This is cheaper than running the rotated addition formula with both inputs equal (12M + 1Ma), but it is a *different* instruction stream. Using it reintroduces an add/double distinction that a power or timing trace can see. The trade-off is spelled out in Section 6.

### 4.4 Completeness condition, stated plainly

| Condition | Consequence |
|---|---|
| `a` is not a cube in `F_p` | no rational points at infinity; rotated addition is complete (works for all inputs) |
| `a` is a cube in `F_p` | rational points at infinity exist; rotated addition has exceptional inputs and needs guards |
| `p ≡ 2 (mod 3)` | every element is a cube, so the first row can never hold; use `p ≡ 1 (mod 3)` |

For a complete formula, therefore, choose `p ≡ 1 (mod 3)` and a non-cube `a`. The cheapest non-cube is a small integer; `a = 2` is a non-cube for roughly two thirds of primes `p ≡ 1 (mod 3)` and costs essentially nothing as a constant multiplication.

### 4.5 The chord (polarisation) formula is not a doubling formula

A tempting alternative derives the third intersection point symbolically. Write `Φ` for the symmetric trilinear form of `F`, set

```
μ = Φ(P₁,P₁,P₂) = aX₁²X₂ + Y₁²Y₂ + Z₁²Z₂ − (d/3)(X₁Y₁Z₂ + X₁Y₂Z₁ + X₂Y₁Z₁)
ν = Φ(P₁,P₂,P₂) = aX₁X₂² + Y₁Y₂² + Z₁Z₂² − (d/3)(X₁Y₂Z₂ + X₂Y₁Z₂ + X₂Y₂Z₁)
```

and take `P₁ + P₂ = [μX₂ − νX₁ : μZ₂ − νZ₁ : μY₂ − νY₁]`.

This is correct for distinct points (the appendix confirms agreement with the rotated formula on 3,980 of 3,980 distinct pairs). But at `P₂ = P₁` both `μ` and `ν` reduce to `Φ(P₁,P₁,P₁) = F(P₁) = 0`, so the output is `(0, 0, 0)` for *every* doubling (confirmed on 200 of 200 points). It is a chord formula, not a unified one, and it is also more expensive (its natural evaluation is well above 12M). It is mentioned here only because it is sometimes presented as "the unified addition law", which it is not. Use the rotated formula of Section 4.2.

---

## 5. Operation counts and comparison with other models

Representative projective-coordinate costs, per [EFD]. Exact numbers depend on coordinate system, mixed vs. full addition, and curve constants; consult the database for the variant you implement.

| Model | Addition | Doubling | One formula for both? |
|---|---|---|---|
| Twisted Hessian, rotated (this document) | 12M + 1Ma | 12M + 1Ma (via addition) or 6M + 3S + 1Ma (dedicated) | **Yes**, when `a` is a non-cube |
| Twisted Hessian, standard | 12M | — | No |
| Twisted Edwards, extended coords [HWCD08] | 9M + 1Ma (8M for `a = −1`, dedicated) | 4M + 4S + 1Ma | Yes (complete when `a` is a square and `d` is a non-square) |
| Short Weierstrass, Jacobian | ≈ 11M + 5S | ≈ 3M + 5S (`a = −3`) | No |
| Short Weierstrass, complete [RCB16] | ≈ 12M + 2Ma (`a = −3`) | ≈ 8M + 3S | Yes |
| Montgomery ladder, x-only | 5M + 4S + 1Ma per ladder step (one add + one double) | | Yes (differential) |

**Honest summary.** Twisted Edwards curves are faster in every column and are also complete. The reasons to pick the twisted Hessian model are narrower:

- The curve you must support already has a rational point of order 3 (e.g. a legacy or standards-mandated curve with cofactor 3), and twisted Edwards form requires a rational point of order 4 instead.
- A single 12M formula with no squarings and no `d` constant simplifies a small fixed-function datapath: one multiplier schedule, no special-casing.
- The addition formula's symmetry (Section 7) gives a clean argument that addition and doubling are indistinguishable at the instruction level.

If none of those apply, implement twisted Edwards.

---

## 6. Scalar multiplication

### 6.1 Uniform ladder

With a complete addition formula, the Montgomery-style ladder needs no case distinctions:

```
Input:  scalar k (n bits, MSB first), point P
R0 ← O ;  R1 ← P
for i from n−1 down to 0:
    b   ← bit i of k
    cswap(R0, R1, b)            # constant-time conditional swap
    R1  ← ADD(R0, R1)           # rotated formula
    R0  ← ADD(R0, R0)           # rotated formula again (doubling)
    cswap(R0, R1, b)
return R0                       # = [k]P
```

Every iteration executes exactly two rotated additions and two swaps, regardless of `k`. Cost per scalar bit: **24M + 2Ma**. For a 256-bit scalar, about 6,200M.

### 6.2 Cheaper variant with dedicated doubling

Replacing `ADD(R0, R0)` by the dedicated doubling formula gives **18M + 3S + 2Ma** per bit (≈ 5,400M for 256 bits), a saving of roughly 20%. The doubling then has a different operation sequence from the addition. Whether an attacker can exploit that depends on what else is in the trace; in a ladder *both* operations run every iteration in the same order, so the operation-type distinction does not by itself reveal the bit. The bit is protected by the swap. The uniform variant is still the more conservative choice, and the one this document recommends when the 20% does not matter.

### 6.3 Fixed-window methods

Window methods (`w = 4` or `5`) reduce the number of additions but require a precomputed table indexed by secret scalar digits. In a constant-time implementation the table must be read with a full linear scan and masked selects, not an indexed load. Window methods also use the dedicated doubling formula between additions, so they inherit the add/double distinction of 6.2. They are appropriate when throughput dominates and the implementation is reviewed for table-access leakage.

---

## 7. Side-channel properties

### 7.1 What the complete formula gives

A scalar multiplication built on one complete addition formula has these properties by construction:

1. **No data-dependent branches.** The formula is a polynomial map; nothing in it tests a field element for zero or compares points.
2. **No exceptional-input branch.** Because the formula never produces `(0, 0, 0)` for valid inputs (Section 4.4), there is no "if degenerate then handle specially" path for an attacker to trigger with a chosen input point.
3. **Identical operation sequence for addition and doubling.** Both are the same 12M + 1Ma schedule reading the same register slots. An attacker who can classify "add vs. double" from a trace gains nothing, because the classifier has no feature to find.

Property 3 is the one usually quoted. Properties 1 and 2 are at least as important in practice, because exceptional-point attacks on incomplete formulas are the classical way to extract a scalar from a Weierstrass implementation.

### 7.2 What it does not give

The complete formula is one layer. It does nothing about:

- **The conditional swap.** `cswap` must be implemented with arithmetic masks, never with a branch or an indexed move. A branching swap leaks every bit of `k` directly.
- **Operand-dependent leakage inside the multiplier.** The intermediate values of the formula depend on the coordinates of `R0` and `R1`, which depend on `k`. A first-order DPA that correlates Hamming weights of intermediates against hypotheses on `k` is **not** prevented by the formula. The standard countermeasure is projective randomisation: multiply `[X : Y : Z]` by a fresh random `r ∈ F_p*` at the start of each scalar multiplication (3M) and optionally between ladder steps. Scalar blinding (`k' = k + m·#E` for random `m`) is a second, independent layer.
- **Timing of the field multiplier.** If the underlying modular multiplication is not constant time, nothing above it is.
- **Memory access patterns.** Register addresses in the formula are fixed, which is good. Precomputation tables (Section 6.3) are not.
- **Fault attacks.** A fault that pushes a point off the curve is not detected by the formula. Validate the input point (`aX³ + Y³ + Z³ − dXYZ = 0`) before and the output point after each scalar multiplication. Cost: a few multiplications, negligible against the ladder.

### 7.3 On "mutual information ≈ 0"

It is sometimes said that a unified formula makes the mutual information between operation type and power trace zero. Read literally, that is a design *goal*, not a theorem. The operation-type channel is removed by construction; other channels remain (Section 7.2). The correct way to make a statement about leakage is to measure it: run a Test Vector Leakage Assessment [GJJR11] on the physical implementation with and without randomisation, and report the t-statistics. Formulas alone do not certify a device.

---

## 8. Parameter requirements

A production curve `TH(a,d)` over `F_p` must satisfy all of the following.

| Requirement | Reason |
|---|---|
| `p` prime, `p > 3`, `p ≡ 1 (mod 3)` | needed so that a non-cube `a` exists |
| `a` not a cube in `F_p` | completeness of the rotated addition formula |
| `a(27a − d³) ≠ 0` | smoothness |
| `#E = h · r` with `r` prime, `h = 3` | **the cofactor is at least 3 for every twisted Hessian curve** (Section 3.5); `h = 3` is the best achievable. `h = 1` is impossible in this model |
| `r ≥ 2^(2λ)` | ρ-method security level `λ` (≈ 256-bit `r` for 128-bit security) |
| `r ∤ p^k − 1` for small `k` | embedding degree large; rules out MOV/Frey–Rück transfer |
| `r ≠ p` | rules out anomalous-curve attack |
| quadratic twist also has a large prime factor | twist security, in case an invalid-curve or fault pushes computation onto the twist |
| `|disc(End(E))|` large | avoid small-CM-discriminant curves unless chosen deliberately |

Point counting is done with the Schoof–Elkies–Atkin algorithm, available in Sage, PARI/GP and Magma; it is polynomial time and takes seconds to minutes at 256 bits. Expect to sample on the order of `log p` candidate `(a, d)` pairs before one has `#E = 3r` with `r` prime and a twist-secure complement.

**Converting an existing curve.** Any Weierstrass curve `y² = x³ + Ax + B` over `F_p` with a rational point `T` of order 3 can be brought to twisted Hessian form over `F_p`; the construction in [BCKL15, Section 5] uses the tangent line at `T`. Conversely, a curve with no rational 3-torsion (for example, every curve of prime order such as the NIST P-curves, Curve25519's prime-order subgroup aside) has no twisted Hessian model over `F_p`. In particular, **the NIST P-256 prime with `a = 1`, `h = 1` is not a valid twisted Hessian parameter set**; the prime is usable, but the group order will be divisible by 3 and the cofactor cannot be 1.

---

## 9. Implementation notes for an HSM

### 9.1 Datapath

- One modular multiplier (Montgomery or Barrett), constant latency, no early termination.
- Register file of at least 14 field elements: 6 inputs, 6 cross terms, 2 accumulators. The outputs can overwrite the cross terms.
- The addition formula's dependency graph has depth 2 (cross terms, then products). With two multipliers the 12M schedule runs in 6 cycles of multiplier latency; with four, in 3. Beyond four, the `6 → 6` structure gives diminishing returns.
- No field inversion inside the ladder. Convert to affine once at the end with a single constant-time inversion (Fermat exponentiation `z^(p−2)`, or a constant-time binary GCD).

### 9.2 Checklist

```
[ ] p ≡ 1 (mod 3); a is a non-cube; a(27a − d³) ≠ 0           (Section 8)
[ ] #E = 3r, r prime, twist secure                           (Section 8)
[ ] rotated addition used for every ladder step              (Section 6.1)
[ ] cswap implemented with masks                             (Section 7.2)
[ ] projective randomisation at the start of each scalar mul (Section 7.2)
[ ] scalar blinding (optional second layer)                  (Section 7.2)
[ ] input and output points validated on the curve           (Section 7.2)
[ ] field multiplier timing independent of operands          (Section 7.2)
[ ] no secret-indexed memory access anywhere                 (Section 6.3)
[ ] TVLA measured on silicon, with and without randomisation (Section 7.3)
[ ] formulas regression-tested against the appendix script   (Section 11)
```

### 9.3 Interoperability

The group is the same abstract group as the Weierstrass curve it came from, so points can be exported in Weierstrass coordinates for protocols that expect them (ECDSA, ECDH with standard encodings). Conversion costs a handful of field operations plus one inversion each way. Keep the twisted Hessian representation internal to the scalar-multiplication core.

---

## 10. Known pitfalls

These are errors that appear in informal descriptions of this model. Each is contradicted by the computation in the appendix.

1. **"The addition law is 12M + 6S."** The rotated addition is 12M + 1Ma with no squarings. The 6S figure comes from an unoptimised evaluation that squares each coordinate separately instead of pairing cross terms. A 17M count arises the same way.
2. **"The polarisation formula is unified."** It returns `(0, 0, 0)` on every doubling (Section 4.5).
3. **"Cofactor 1 is achievable."** Every twisted Hessian curve has order divisible by 3 (Section 3.5).
4. **"`j = d³(d³ − 216a)³ / …`."** The sign is `+216a` (Section 3.2).
5. **"The unified formula proves I(op; trace) = 0."** It removes one channel; the rest must be measured (Section 7.3).
6. **"The group law can be computed with 16-bit fixed-point shift-and-add arithmetic."** The group law is modular arithmetic over a 256-bit prime field. Fixed-point approximations are not applicable; every operation must be an exact modular multiplication.
7. **"Twisted Hessian is faster than twisted Edwards."** It is not (Section 5). Its advantages are structural, not speed.

---

## 11. Verification

The script in the appendix runs in a few seconds with no dependencies beyond Python 3. Over `p = 1009` (`p ≡ 1 mod 3`), `a = 2` (a non-cube), `d = 5`, it:

- enumerates all 975 rational points and checks `#E ≡ 0 (mod 3)` and the Hasse bound;
- checks that the rotated addition formula sends every sampled pair of points to a point on the curve, satisfies `P + O = P`, `P + (−P) = O`, and `(P + P) − P = P`;
- checks the rotated formula against the standard formula on all pairs where the standard formula is defined (agreement on every pair);
- checks the dedicated doubling formula against `ADD(P, P)` (agreement on every point);
- checks associativity on 300 random triples;
- finds **zero** pairs out of 40,000 on which the rotated formula returns `(0, 0, 0)`;
- counts rational 3-torsion (3 points: `O` and two others, as predicted for `p ≡ 1 mod 3`);
- shows the polarisation formula agrees on distinct points and fails on every doubling;
- confirms `#E ≡ 0 (mod 3)` for 24 further parameter sets over `p ≡ 2 (mod 3)`, and locates the `j = 0` curves at `d = 0` and `d = −6` by a supersingularity test, fixing the sign in the `j`-invariant.

Running it is the acceptance test for any port of these formulas to another language or to RTL: generate the same point set, apply the ported formula, compare coordinates after normalisation.

---

## 12. Limitations and open items

- **No hardware measurements are reported here.** Operation counts are exact; cycle counts, area and leakage figures depend on the datapath and must be measured.
- **No production parameters are given.** Section 8 lists the constraints; generating and certifying a curve is a separate task with its own review.
- **Cofactor 3 is a protocol consideration.** Protocols that assume prime-order groups (some signature and key-exchange proofs do) need either cofactor clearing (multiply by 3) or a re-analysis.
- **Classical security only.** The discrete-logarithm problem on these curves offers no resistance to Shor's algorithm. This model is a choice of *implementation* for classical ECC, not a post-quantum scheme.
- The connection between this model and anything outside elliptic-curve arithmetic (machine-learning optimisation, renormalisation, information geometry) is not established by anything in this document, and no such connection is claimed.

---

## 13. References

- **[BCKL15]** D. J. Bernstein, C. Chuengsatiansup, D. Kohel, T. Lange. *Twisted Hessian curves.* LATINCRYPT 2015, LNCS 9230, pp. 269–294. Full version: https://eprint.iacr.org/2015/781
- **[BL07]** D. J. Bernstein, T. Lange. *Faster addition and doubling on elliptic curves.* ASIACRYPT 2007, LNCS 4833, pp. 29–50.
- **[EFD]** D. J. Bernstein, T. Lange. *Explicit-Formulas Database.* https://hyperelliptic.org/EFD/
- **[FJ10]** R. R. Farashahi, M. Joye. *Efficient arithmetic on Hessian curves.* PKC 2010, LNCS 6056, pp. 243–260.
- **[GJJR11]** G. Goodwill, B. Jun, J. Jaffe, P. Rohatgi. *A testing methodology for side-channel resistance validation.* NIST Non-Invasive Attack Testing Workshop, 2011.
- **[HWCD08]** H. Hisil, K. K.-H. Wong, G. Carter, E. Dawson. *Twisted Edwards curves revisited.* ASIACRYPT 2008, LNCS 5350, pp. 326–343.
- **[JQ01]** M. Joye, J.-J. Quisquater. *Hessian elliptic curves and side-channel attacks.* CHES 2001, LNCS 2162, pp. 402–410.
- **[KJJ99]** P. Kocher, J. Jaffe, B. Jun. *Differential power analysis.* CRYPTO 1999, LNCS 1666, pp. 388–397.
- **[Mon87]** P. L. Montgomery. *Speeding the Pollard and elliptic curve methods of factorization.* Mathematics of Computation 48 (1987), pp. 243–264.
- **[RCB16]** J. Renes, C. Costello, L. Batina. *Complete addition formulas for prime order elliptic curves.* EUROCRYPT 2016, LNCS 9665, pp. 403–428.
- **[Sil09]** J. H. Silverman. *The Arithmetic of Elliptic Curves*, 2nd ed. GTM 106, Springer, 2009.
- **[Sma01]** N. P. Smart. *The Hessian form of an elliptic curve.* CHES 2001, LNCS 2162, pp. 118–125.

---

## Appendix: verification script

```python
# th_check.py — numeric checks for every formula in this document.
# Python 3, no dependencies. Runs in a few seconds.

import random

p = 1009                     # p ≡ 1 (mod 3): non-cubes exist
assert p % 3 == 1

def is_cube(v):
    return any(pow(t, 3, p) == v % p for t in range(p))

a = next(v for v in range(2, p) if not is_cube(v))   # a = 2 for p = 1009
d = 5
assert (a * (27 * a - d**3)) % p != 0                # smooth

def on_curve(P):
    X, Y, Z = P
    return (a*X**3 + Y**3 + Z**3 - d*X*Y*Z) % p == 0

# All rational points: affine chart plus points at infinity (Z = 0).
pts = [(x, y, 1) for x in range(p) for y in range(p)
       if (a*x**3 + y**3 + 1 - d*x*y) % p == 0]
inf = sum(1 for y in range(1, p) if (a + y**3) % p == 0)   # X = 1, Z = 0
N = len(pts) + inf
print("affine:", len(pts), " at infinity:", inf, " #E =", N)
assert N % 3 == 0, "group order must be divisible by 3"
assert abs(N - (p + 1)) <= 2 * int(p**0.5) + 1, "Hasse bound"
assert inf == 0, "non-cube a => no rational points at infinity"

def norm(P):
    X, Y, Z = [v % p for v in P]
    for c in (Z, Y, X):
        if c:
            inv = pow(c, p - 2, p)
            return (X*inv % p, Y*inv % p, Z*inv % p)
    raise ValueError("zero point")

O = norm((0, -1, 1))

def neg(P):
    X, Y, Z = P
    return norm((X, Z, Y))

def add_std(P, Q):                       # 12M, incomplete
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    return norm((X1*X1*Y2*Z2 - X2*X2*Y1*Z1,
                 Z1*Z1*X2*Y2 - Z2*Z2*X1*Y1,
                 Y1*Y1*X2*Z2 - Y2*Y2*X1*Z1))

def add_rot(P, Q):                       # 12M + 1Ma, complete for non-cube a
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    return norm((Z2*Z2*X1*Z1 - Y1*Y1*X2*Y2,
                 Y2*Y2*Y1*Z1 - a*X1*X1*X2*Z2,
                 a*X2*X2*X1*Y1 - Z1*Z1*Y2*Z2))

def dbl(P):                              # 6M + 3S + 1Ma
    X1, Y1, Z1 = P
    return norm((X1*(Z1**3 - Y1**3),
                 Z1*(Y1**3 - a*X1**3),
                 Y1*(a*X1**3 - Z1**3)))

random.seed(1)
S = [norm(P) for P in random.sample(pts, 200)]

# 1. Standard formula: defined on most distinct pairs, never on doublings.
std_fail = 0
for P in S:
    for Q in S[:30]:
        try:
            assert on_curve(add_std(P, Q))
        except ValueError:
            std_fail += 1
print("standard-law undefined pairs (expected > 0):", std_fail)

# 2. Rotated formula: group axioms, agreement with standard law, doubling.
for P in S:
    for Q in S[:20]:
        if P != Q:
            try:
                Rs = add_std(P, Q)
            except ValueError:
                continue
            assert add_rot(P, Q) == Rs
    assert add_rot(P, O) == P and add_rot(O, P) == P
    assert add_rot(P, neg(P)) == O
    D = add_rot(P, P)
    assert on_curve(D) and D == dbl(P)
    assert add_rot(D, neg(P)) == P

# 3. Associativity.
for _ in range(300):
    P, Q, R = random.sample(S, 3)
    assert add_rot(add_rot(P, Q), R) == add_rot(P, add_rot(Q, R))

# 4. Completeness: rotated formula never returns (0:0:0).
bad = 0
for P in S:
    for Q in S:
        X1, Y1, Z1 = P; X2, Y2, Z2 = Q
        out = ((Z2*Z2*X1*Z1 - Y1*Y1*X2*Y2) % p,
               (Y2*Y2*Y1*Z1 - a*X1*X1*X2*Z2) % p,
               (a*X2*X2*X1*Y1 - Z1*Z1*Y2*Z2) % p)
        if out == (0, 0, 0):
            bad += 1
print("rotated-law exceptional pairs (expected 0):", bad)
assert bad == 0

# 5. Rational 3-torsion.
tors = sum(1 for P in (norm(q) for q in pts) if add_rot(dbl(P), P) == O)
print("rational 3-torsion points incl. O (expected 3 for p ≡ 1 mod 3):", tors)

# 6. The chord / polarisation formula is NOT a doubling formula.
inv3 = pow(3, p - 2, p)
def add_munu(P, Q):
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    mu = (a*X1*X1*X2 + Y1*Y1*Y2 + Z1*Z1*Z2
          - d*inv3*(X1*Y1*Z2 + X1*Y2*Z1 + X2*Y1*Z1))
    nu = (a*X1*X2*X2 + Y1*Y2*Y2 + Z1*Z2*Z2
          - d*inv3*(X1*Y2*Z2 + X2*Y1*Z2 + X2*Y2*Z1))
    return ((mu*X2 - nu*X1) % p, (mu*Z2 - nu*Z1) % p, (mu*Y2 - nu*Y1) % p)
agree = sum(1 for P in S for Q in S[:20]
            if P != Q and norm(add_munu(P, Q)) == add_rot(P, Q))
total = sum(1 for P in S for Q in S[:20] if P != Q)
zeros = sum(1 for P in S if add_munu(P, P) == (0, 0, 0))
print(f"chord formula: agrees on {agree}/{total} distinct pairs,"
      f" returns (0:0:0) on {zeros}/{len(S)} doublings")

# 7. Group order ≡ 0 (mod 3) for other parameters, and j = 0 location.
#    Over p ≡ 2 (mod 3), j = 0 curves are supersingular: #E = p + 1.
for q in (1013, 1019):
    assert q % 3 == 2
    for (aa, dd) in ((1, 5), (2, 3), (7, 11), (1, 0), (1, -6), (1, 6)):
        if (aa * (27*aa - dd**3)) % q == 0:
            continue
        n = sum(1 for x in range(q) for y in range(q)
                if (aa*x**3 + y**3 + 1 - dd*x*y) % q == 0)
        n += sum(1 for y in range(1, q) if (aa + y**3) % q == 0)
        assert n % 3 == 0
        if (aa, dd) in ((1, 0), (1, -6)):
            assert n == q + 1, "expected j = 0 (supersingular)"
        if (aa, dd) == (1, 6):
            assert n != q + 1, "d = +6 is not j = 0"
print("order ≡ 0 mod 3 and j = 0 at d ∈ {0, −6} (a = 1): confirmed")
print("ALL CHECKS PASSED   p =", p, " a =", a, " d =", d)
```

Expected output:

```
affine: 975  at infinity: 0  #E = 975
standard-law undefined pairs (expected > 0): 43
rotated-law exceptional pairs (expected 0): 0
rational 3-torsion points incl. O (expected 3 for p ≡ 1 mod 3): 3
chord formula: agrees on 3980/3980 distinct pairs, returns (0:0:0) on 200/200 doublings
order ≡ 0 mod 3 and j = 0 at d ∈ {0, −6} (a = 1): confirmed
ALL CHECKS PASSED   p = 1009  a = 2  d = 5
```
