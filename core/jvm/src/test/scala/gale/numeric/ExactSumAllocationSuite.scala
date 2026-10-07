package gale.numeric

import java.lang.management.ManagementFactory

/** Engineering allocation probe on supported HotSpot; not a workload/RSS benchmark. */
class ExactSumAllocationSuite extends munit.FunSuite:
  test("qualified JVM ratio allocation fits the conservative payload estimate") {
    val host = ManagementFactory.getThreadMXBean
    assume(host.isInstanceOf[com.sun.management.ThreadMXBean], "thread allocation measurement is unavailable")
    val bean = host.asInstanceOf[com.sun.management.ThreadMXBean]
    assume(bean.isThreadAllocatedMemorySupported, "thread allocation measurement is unavailable")
    val prior = bean.isThreadAllocatedMemoryEnabled
    if !prior then bean.setThreadAllocatedMemoryEnabled(true)
    try {
      val cases = Vector(
        (Double.MaxValue, 1.0),
        (Double.MinPositiveValue, Double.MaxValue),
        (Double.MaxValue, Double.MinPositiveValue),
        (-Double.MaxValue, -1.0),
        (1.0, Double.MaxValue),
        (Double.MinPositiveValue, 2.0)
      ).map { (n, d) =>
        val a = ExactSum.zero(); val b = ExactSum.zero()
        assert(a.add(n).isRight); assert(b.add(d).isRight)
        (0 until 60).foreach { _ => assert(a.addAll(a).isRight); assert(b.addAll(b).isRight) }
        (a, b)
      }
      cases.foreach { (a, b) => (0 until 10).foreach(_ => a.ratio(b)) }
      val id = Thread.currentThread().threadId()
      cases.foreach { (a, b) =>
        val before = bean.getThreadAllocatedBytes(id)
        val result = a.ratio(b)
        val allocated = bean.getThreadAllocatedBytes(id) - before
        assert(result.isRight)
        assert(allocated >= 0)
        assert(allocated <= a.resources.ratioScratchBytes, s"$allocated > ${a.resources.ratioScratchBytes}")
      }
    } finally {
      if !prior then bean.setThreadAllocatedMemoryEnabled(false)
    }
  }
