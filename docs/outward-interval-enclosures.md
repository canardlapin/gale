# Finite outward interval enclosures

`gale.linalg.OutwardInterval` is a small portable certification primitive. Its
closed finite endpoints are constructed with `Either`; it rejects NaN,
infinity, inverted bounds, a divisor containing zero, and finite operations
that overflow. Addition, subtraction, multiplication and division expand each
ordinary IEEE-754 endpoint result outwards with `Math.nextDown` and
`Math.nextUp`. `square`, `streamingDot`, and `streamingAbsoluteSum` use those
same operations and retain only scalar state while streaming arrays.

The scope is intentional. This is not an eigenvalue, factorization, inverse,
matrix, or HRF certificate. It has no transcendental operations: in particular
it does not expose square root because a portable correctly-rounded/sound
outward expansion proof has not been supplied. Consumers requiring family,
whitening, reconstruction, or spectral envelopes must bind separately proved
enclosures and refuse when those capabilities are unavailable.

The focused controls use an exact dyadic-`BigInt` oracle for nondegenerate
multiplication rectangles in each sign quadrant and exact rational corners for
strictly positive and strictly negative divisors. They also cover one-sided
squares, refusal when expanding an endpoint adjacent to `Double.MaxValue`,
empty streams, late NaN/infinity input, and signed-zero endpoints. The stream
implementation has constant live scalar state by source inspection; these
tests do not measure allocations or claim that allocations are zero.

The result encloses the real result of its finite primitive inputs under the
ordinary binary64 model. It does not turn a host runtime, an input source, or a
scientific workflow into a qualified certificate.
