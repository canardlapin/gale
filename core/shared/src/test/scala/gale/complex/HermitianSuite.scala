package gale.complex

import munit.FunSuite

class HermitianSuite extends FunSuite:
  private val budget=BigInt(64)*1024*1024
  private val policy=HermitianPolicy(HermitianCutoff.Relative(1e-10))
  private def ok[A](v: Either[?,A]): A=v.fold(e=>fail(e.toString),identity)
  private def matrix(n: Int,m: Int,r: Vector[Double],i: Vector[Double]): OwnedComplex=ok(OwnedComplex.copying(n,m,r.toArray,i.toArray,budget))
  private def compare(a: OwnedComplex,r: Vector[Double],i: Vector[Double],tol: Double=1e-9): Unit =
    a.realCopy.toVector.zip(r).foreach((v,e)=>assertEqualsDouble(v,e,tol))
    a.imaginaryCopy.toVector.zip(i).foreach((v,e)=>assertEqualsDouble(v,e,tol))
  for (f,index)<-HermitianReference.cases.zipWithIndex do
    test(s"independent NumPy Hermitian eigensystem factor projector and solve $index"):
      val a=matrix(f.n,f.n,f.real,f.imaginary);val result=ok(Hermitian.factor(a,policy,budget))
      assertEquals(result.record.rank,f.rank)
      result.record.eigenvalues.zip(f.eigenvalues).foreach((v,e)=>assertEqualsDouble(v,e,1e-10))
      assert(result.certificate.eigenResidual<=1e-9);assert(result.certificate.orthogonalityError<=1e-9)
      compare(ok(result.projector(budget)),f.projectorReal,f.projectorImaginary)
      compare(ok(result.pseudoInverse(budget)),f.inverseReal,f.inverseImaginary)
      val factor=ok(result.psdFactor(budget));val reconstructed=ok(factor.multiply(ok(factor.conjugateTranspose(budget)),budget))
      compare(reconstructed,f.real,f.imaginary,1e-9)
      val b=matrix(f.n,2,f.rhsReal,f.rhsImaginary)
      val solved=ok(result.solve(b,HermitianSolveMode.RetainedLeastSquares,budget))
      compare(solved.values,f.solutionReal,f.solutionImaginary)
      if f.rank==f.n then
        assert(ok(result.solve(b,HermitianSolveMode.FullRank,budget)).relativeResidual<=1e-9)
      else
        assertEquals(result.solve(b,HermitianSolveMode.FullRank,budget),Left(ComplexFailure.RankDeficient(f.rank,f.n)))
        assert(solved.relativeDiscardedRhs>0.01)

  test("conjugate transpose product and all public copies preserve immutable ownership"):
    val re=Array(1.0,2.0,3.0,4.0);val im=Array(0.0,1.0,-2.0,0.0)
    val a=ok(OwnedComplex.copying(2,2,re,im,budget));re(0)=99;im(1)=99
    val t=ok(a.conjugateTranspose(budget));compare(t,Vector(1,3,2,4),Vector(0,2,-1,0))
    val copy=a.realCopy;copy(0)=99;assertEquals(a(0,0).re,1.0)
    val gram=ok(a.multiply(t,budget));compare(gram,Vector(6,11,11,29),Vector(0,6,-6,0))
    assert(a.multiply(matrix(1,1,Vector(1),Vector(0)),budget).isLeft)
    assert(a.conjugateTranspose(a.storageBytes-1).isLeft)
    assert(a.multiply(t,a.storageBytes-1).isLeft)

  test("zero repeated eigenvalues and explicit absolute cutoff retain declared rank"):
    val zero=matrix(3,3,Vector.fill(9)(0),Vector.fill(9)(0))
    val z=ok(Hermitian.factor(zero,policy,budget));assertEquals(z.record.rank,0)
    compare(ok(z.pseudoInverse(budget)),Vector.fill(9)(0),Vector.fill(9)(0))
    val identity=matrix(3,3,Vector(1,0,0,0,1,0,0,0,1),Vector.fill(9)(0))
    assertEquals(ok(Hermitian.factor(identity,policy,budget)).record.rank,3)
    val cutoff=policy.copy(cutoff=HermitianCutoff.Absolute(1.0))
    assertEquals(ok(Hermitian.factor(identity,cutoff,budget)).record.rank,0)
    val diag=matrix(2,2,Vector(-1e-14,0,0,1),Vector.fill(4)(0))
    val accepted=ok(Hermitian.factor(diag,policy,budget))
    assertEquals(accepted.record.acceptedNegativeRoundoff,Vector(-1e-14))
    assert(Hermitian.factor(diag,policy.copy(psdTolerance=0),budget).isLeft)

  test("large and small scale preserve eigenvectors rank and rescaled eigenvalues"):
    val f=HermitianReference.cases.head
    for scale<-Vector(1e-150,1e150) do
      val a=matrix(f.n,f.n,f.real.map(_*scale),f.imaginary.map(_*scale))
      val result=ok(Hermitian.factor(a,policy,budget));assertEquals(result.record.rank,f.rank)
      result.record.eigenvalues.zip(f.eigenvalues).foreach((v,e)=>assertEqualsDouble(v/scale,e,1e-10))
      compare(ok(result.projector(budget)),f.projectorReal,f.projectorImaginary)

  test("indefinite non-Hermitian nonfinite and unsupported shape inputs refuse"):
    val bad=Vector(
      matrix(2,2,Vector(1,0,0,1),Vector(0,2,-2,0)),
      matrix(2,2,Vector(1,1,0,1),Vector.fill(4)(0)),
      matrix(2,2,Vector(1,0,0,1),Vector(1,0,0,0)),
      matrix(2,1,Vector(1,1),Vector(0,0)))
    bad.foreach(a=>assert(Hermitian.factor(a,policy,budget).isLeft))
    assert(OwnedComplex.copying(1,1,Array(Double.NaN),Array(0.0),budget).isLeft)
    assert(OwnedComplex.copying(1,1,Array(0.0),Array(Double.PositiveInfinity),budget).isLeft)
    assert(OwnedComplex.copying(1,2,Array(1.0),Array(0.0),budget).isLeft)
    assert(OwnedComplex.inspect(Int.MaxValue,Int.MaxValue,BigInt(Long.MaxValue)*16).isLeft)
    assert(Hermitian.inspect(129,policy,budget).isLeft)

  test("admitted roundoff projection is recorded and strict symmetry can refuse it"):
    val a=matrix(2,2,Vector(2,0,0,2),Vector(0,1+1e-14,-1,0))
    val r=ok(Hermitian.factor(a,policy,budget));assert(r.certificate.relativeHermitianCorrection>0)
    assert(Hermitian.factor(a,policy.copy(symmetryTolerance=0),budget).isLeft)

  test("budget boundaries cancellation and exhausted sweeps never publish a factor"):
    val f=HermitianReference.cases(3);val a=matrix(f.n,f.n,f.real,f.imaginary)
    val bytes=Hermitian.numericalBytes(f.n)
    assertEquals(Hermitian.factor(a,policy,bytes-1),Left(ComplexFailure.Budget(bytes,bytes-1)))
    assert(Hermitian.factor(a,policy,bytes).isRight)
    assertEquals(Hermitian.factor(a,policy,budget,()=>true),Left(ComplexFailure.Cancelled))
    var calls=0
    assertEquals(Hermitian.factor(a,policy,budget,()=>{calls+=1;calls>5}),Left(ComplexFailure.Cancelled))
    assert(Hermitian.factor(a,policy.copy(maxSweeps=1),budget) match
      case Left(ComplexFailure.NonConverged(1,_))=>true
      case _=>false)
    for p<-Vector(policy.copy(maxSweeps=0),policy.copy(cutoff=HermitianCutoff.Relative(Double.NaN)),
      policy.copy(convergenceTolerance=0),policy.copy(certificateTolerance=1)) do
      assert(Hermitian.factor(a,p,budget).isLeft)
