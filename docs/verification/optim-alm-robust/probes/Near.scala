import gale.optim.*
import gale.linalg.{DMat,DVec}
object Near:
 def main(args:Array[String]):Unit =
  def pt(x:Double,y:Double)=DMat.dense(2,1,Seq(x,y))
  val c=NonlinearConstraints(2,2,0)(x=>Right(ConstraintEvaluation(DVec.fromSeq(Seq(x(0,0)+.001*x(1,0)-1,x(0,0)-.001*x(1,0)-1)),DMat.dense(2,2,Seq(1,.001,1,-.001))))).toOption.get
  val f=DifferentiableObjective(2)(x=>{val dx=x(0,0)-3;val dy=x(1,0)-2;Right(ObjectiveEvaluation(.5*(dx*dx+dy*dy),pt(dx,dy)))}).toOption.get
  for adapt<-Seq(false,true);reuse<-Seq(false,true) do
   val r=AugmentedLagrangian.minimize(f,c,pt(0,0),config=AugmentedLagrangianConfig(adaptiveInnerTolerance=adapt,reuseCurvatureScale=reuse,control=SolverControl(traceCapacity=100)))
   println(s"adaptive=$adapt reuse=$reuse result="+r.map(z=>(z.status,z.primal(0,0),z.primal(1,0),z.evaluations.callbacks,z.diagnostics,z.innerTrace.map(t=>(t.penalty,t.initialCurvatureScale,t.iterations,t.status)))))
