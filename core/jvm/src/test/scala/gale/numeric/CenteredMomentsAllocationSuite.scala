package gale.numeric

class CenteredMomentsAllocationSuite extends munit.FunSuite:
  private def ok[E, A](value: Either[E, A]): A = value.fold(e => fail(e.toString), identity)
  test("ordinary eight-channel exact observations avoid repeated padded-integer limb allocation"):
    val rows = Vector.tabulate(512)(i => Array.tabulate(8)(r => (r + 1) * 1e-5 + i * 1e-10))
    val warm = ok(CenteredMoments.make(8))
    for _ <- 0 until 3; row <- rows do ok(warm.add(row))
    val state = ok(CenteredMoments.make(8))
    val bean = java.lang.management.ManagementFactory.getThreadMXBean.asInstanceOf[com.sun.management.ThreadMXBean]
    assert(bean.isThreadAllocatedMemorySupported)
    bean.setThreadAllocatedMemoryEnabled(true)
    val thread = Thread.currentThread().threadId()
    val before = bean.getThreadAllocatedBytes(thread)
    rows.foreach(row => ok(state.add(row)))
    val allocated = bean.getThreadAllocatedBytes(thread) - before
    assert(allocated < 9437184L, s"512 exact8-channel observations allocated $allocated bytes")
    assertEquals(state.frequency, 512L)
    assert(ok(state.finish()).mean(0).isFinite)
