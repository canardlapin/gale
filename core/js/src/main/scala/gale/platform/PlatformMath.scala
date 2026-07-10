package gale.platform

/** Platform arithmetic primitives used by shared hot loops. */
private[gale] object PlatformMath:
  /** JavaScript has no portable scalar FMA; keep the shared call shape while
    * lowering to the ordinary two-operation expression.
    */
  inline def fma(a: Double, b: Double, c: Double): Double =
    a * b + c
