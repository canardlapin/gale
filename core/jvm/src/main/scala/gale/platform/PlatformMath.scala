package gale.platform

/** Platform arithmetic primitives used by shared hot loops. */
private[gale] object PlatformMath:
  /** One-rounding fused multiply-add on the JVM. */
  inline def fma(a: Double, b: Double, c: Double): Double =
    java.lang.Math.fma(a, b, c)
