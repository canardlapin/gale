package gale.bench
import gale.linalg.{DMat, DVec}
import gale.optim.*
object AllocationProbe:
  def main(args: Array[String]): Unit =
    val n = 128
    val point = DMat.zeros(n, 1)
    val jac = DMat.zeros(1, n)
    val f = DifferentiableObjective(n)(_ => Right(ObjectiveEvaluation(0.0, point))).toOption.get
    val c = NonlinearConstraints(n, 1, 0)(_ => Right(ConstraintEvaluation(DVec.zeros(1), jac))).toOption.get
    def run(): Unit =
      val r = AugmentedLagrangian.minimize(f, c, point).toOption.get
      require(r.status == AugmentedLagrangianStatus.Converged && r.evaluations.callbacks == 2)
    for _ <- 0 until 3000 do run()
    val bean = java.lang.management.ManagementFactory.getThreadMXBean.asInstanceOf[com.sun.management.ThreadMXBean]
    require(bean.isThreadAllocatedMemorySupported)
    bean.setThreadAllocatedMemoryEnabled(true)
    val id = Thread.currentThread().threadId()
    val before = bean.getThreadAllocatedBytes(id)
    for _ <- 0 until 1000 do run()
    val bytes = bean.getThreadAllocatedBytes(id) - before
    println(s"ALLOCATION bytes_per_zero_solve=${bytes / 1000.0} callbacks=2 n=$n")
