package gale.numeric

class CenteredMomentsFixedAllocationSuite extends munit.FunSuite:
  test("4096 admitted exact weighted updates allocate no observation-sized limbs or arrays"):
    val bean=java.lang.management.ManagementFactory.getThreadMXBean.asInstanceOf[com.sun.management.ThreadMXBean]
    assert(bean.isThreadAllocatedMemorySupported)
    if !bean.isThreadAllocatedMemoryEnabled then bean.setThreadAllocatedMemoryEnabled(true)
    val input=Array(1e-6,-2e-6)
    val warm=CenteredMoments.make(2).toOption.get
    var i=0
    while i<20000 do {val _ = warm.add(input,1+(i%7));i+=1}
    val state=CenteredMoments.make(2).toOption.get
    val thread=Thread.currentThread().threadId()
    val before=bean.getThreadAllocatedBytes(thread)
    i=0
    while i<4096 do {val _ = state.add(input,1+(i%7));i+=1}
    val bytes=bean.getThreadAllocatedBytes(thread)-before
    assert(bytes<65536,s"exact updates allocated $bytes bytes; fixed retained limbs required")
    assert(state.finish(0).isRight)
