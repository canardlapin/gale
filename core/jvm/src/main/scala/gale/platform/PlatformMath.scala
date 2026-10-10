package gale.platform

/** Platform arithmetic primitives used by shared hot loops. */
private[gale] object PlatformMath:
  /** One-rounding fused multiply-add on the JVM. */
  inline def fma(a: Double, b: Double, c: Double): Double =
    java.lang.Math.fma(a, b, c)

  /** True where an in-order compare scan beats a `math.max`/`math.min`
    * reduction for extreme values. HotSpot C2 on x86 lowers double
    * `Math.max`/`Math.min` to a compare-and-branch ladder that takes a branch
    * per element and does not vectorize, while a `v <= m` scan with the update
    * out of line takes none on the common path. On AArch64 they are single
    * `fmax`/`fmin` instructions and the reduction is faster.
    */
  val scanExtremes: Boolean =
    System.getProperty("os.arch", "") match
      case "amd64" | "x86_64" | "x86" | "i386" | "i686" => true
      case _                                             => false
