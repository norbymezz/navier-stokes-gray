# Discrete structure: changes, accumulation, scale and representation

This note records a research track that emerged from the Gray-code experiments.
It is deliberately separated from the fluid equations: the objective is to identify reusable discrete structure first, then test whether it reduces cost or clarifies scale coupling in acoustics and fluid mechanics.

## 1. Observable state is not necessarily complete state

Two paths can reach the same displayed coordinate while retaining different orientation or history. We therefore distinguish an observable coordinate x_n from a richer state

S_n = (x_n, g_n, chi_n, W_n)

where g_n is a Gray state, chi_n is a parity/orientation class, and W_n is an operator word recording the path.

Thus x_A = x_B does not imply S_A = S_B.

## 2. Local difference and accumulation

For any sequence a_n, define its cumulative sum

S_n = sum_{k=1..n} a_k.

Then the finite difference recovers the added term:

Delta S_n = S_n - S_(n-1) = a_n.

This gives a discrete structural pair:

local change <-> accumulation
Delta <-> Sigma

For power sums,

S_p(n) = sum_{k=1..n} k^p

we have

S_1(n) = n(n+1)/2

S_2(n) = n(n+1)(2n+1)/6 = S_1(n)(2n+1)/3

S_3(n) = [S_1(n)]^2

So these levels are not independent. In particular, when n = 3m + 1,

S_2(n) / S_1(n) = 2m + 1

and the scale factor is an odd integer.

The general power-sum family belongs to the Faulhaber/Bernoulli framework. This is a reference point for later work; no new identity is claimed here.

## 3. Gray transform as local XOR change

For a fixed-width non-negative binary word b, reflected binary Gray code is

g = b XOR (b >> 1).

In bits,

g_(m-1) = b_(m-1)
g_i = b_(i+1) XOR b_i

This can be read as a local XOR difference along adjacent bit positions.

The inverse is cumulative XOR:

b_(m-1) = g_(m-1)
b_i = b_(i+1) XOR g_i

The analogy with Delta <-> Sigma is structural: Gray stores local change information and the inverse reconstructs by accumulation. It is not an identification of XOR calculus with ordinary real-valued calculus.

## 4. Which bit changes: the ruler sequence

Successive reflected Gray states differ in exactly one bit. If bit positions are numbered from 1 at the least-significant bit, the changed-bit sequence is

1, 2, 1, 3, 1, 2, 1, 4, ...

For transition number k >= 1,

r(k) = nu_2(k) + 1

where nu_2(k) is the exponent of 2 dividing k.

Finite prefixes satisfy the recursive reflection rule

D_n = D_(n-1), n, D_(n-1).

This is the dyadic/ruler structure behind the reflected Gray traversal.

## 5. Half + symmetry -> complete structure

For triangular numbers,

2 S_1(n) = n(n+1).

The triangular half is completed by the pairing

k <-> n + 1 - k.

For reflected Gray code,

G_n = 0 G_(n-1) || 1 reverse(G_(n-1)).

The second half is generated from the first by reflection/reversal plus one new leading bit.

Research question: when a computation contains a comparable half + symmetry/operator structure, can one compute or store only the minimal half plus the reconstruction rule?

## 6. Representation can hide structure

Looking along the cube diagonal n = (1,1,1)/sqrt(3) projects the three mutually orthogonal coordinate axes to three planar directions separated by 120 degrees.

Standing warning for the project:

A property of a representation is not automatically a property of the underlying object.

## 7. Relative-turn path experiment

A path is specified by lengths and relative turns. If the current orientation is theta_j and the next relative turn is delta_j,

theta_(j+1) = theta_j + delta_j

p_(j+1) = p_j + L_(j+1) (cos(theta_(j+1)), sin(theta_(j+1))).

The current exploratory family uses descending lengths such as

8, 7, 6, 5, 4, 3, 2, 1

and relative turns chosen from 45 and 90 degrees, with orientation/sign treated as separate information. Closure is to be tested numerically rather than assumed.

## 8. What to search next

The next mathematical comparison should focus on the intersection of:

- reflected Gray code, ruler sequences and 2-adic/dyadic structure;
- finite differences, power sums, Faulhaber polynomials and Bernoulli numbers;
- Haar/Walsh representations and dyadic multiresolution;
- reversible operator words and minimal state/history information.

The later physical question is narrower: can any of these structures provide a useful coarse-graining or reconstruction map between scales in acoustics or fluid mechanics while preserving the observables required by the benchmark?

Any useful result must be measured as either lower error, lower computational cost for comparable error, lower storage, or a clearer reversible mapping.
