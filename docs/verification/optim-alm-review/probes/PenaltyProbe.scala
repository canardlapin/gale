import gale.optim.*
import gale.linalg.{DMat,DVec}
object Probe:
 def main(args:Array[String]):Unit =
  def pt(x:Double)=DMat.dense(1,1,Seq(x))
  val c=NonlinearConstraints(1,1,0)(x=>Right(ConstraintEvaluation(DVec.fromSeq(Seq(x(0,0))),pt(1)))).toOption.get
  for rho <- Seq(10.0,1e3,1e6,1e9) do
   var calls=0
   val f=DifferentiableObjective(1)(x=>{calls+=1;val d=x(0,0)-1;Right(ObjectiveEvaluation(0.5*d*d,pt(d)))}).toOption.get
   val r=AugmentedLagrangian.minimize(f,c,pt(0),config=AugmentedLagrangianConfig(maxIterations=1,initialPenalty=rho,maximumPenalty=1e12)).toOption.get
   println(s"rho=$rho calls=$calls inner=${r.innerIterations} status=${r.status} x=${r.primal(0,0)} stationarity=${r.diagnostics.stationarity}")
