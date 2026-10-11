package gale.numeric

class CenteredMomentsFixedLimbsSuite extends munit.FunSuite:
  private def ok[A](value: Either[CenteredMomentsError,A]): A=value.fold(e=>fail(e.message),identity)
  private def units(value: Double): BigInt =
    val bits=java.lang.Double.doubleToRawLongBits(value)
    val exponent=((bits>>>52)&0x7ffL).toInt
    val fraction=bits&0xfffffffffffffL
    val magnitude=if exponent==0 then BigInt(fraction) else BigInt(fraction|(1L<<52))<<(exponent-1)
    if bits<0 then -magnitude else magnitude
  private def verify(samples: Vector[(Array[Double],Long)]): Unit =
    val n=samples.head._1.length
    val state=ok(CenteredMoments.make(n));samples.foreach((values,weight)=>ok(state.add(values,weight)))
    val mass=samples.map(_._2).sum
    val sums=Array.tabulate(n)(r=>samples.foldLeft(BigInt(0))((sum,p)=>sum+units(p._1(r))*p._2))
    val products=Array.tabulate(n,n)((r,c)=>samples.foldLeft(BigInt(0))((sum,p)=>sum+units(p._1(r))*units(p._1(c))*p._2))
    val expectedMeans=sums.map(s=>ExactDyadic.roundRatio(s,BigInt(mass)<<1074))
    val expected=Array.tabulate(n,n)((r,c)=>ExactDyadic.roundRatio(BigInt(mass)*products(r)(c)-sums(r)*sums(c),BigInt(mass)*mass<<2148))
    val result=state.finish(0)
    if expected.flatten.exists(!_.isFinite) then assert(result.isLeft)
    else
      val summary=ok(result)
      for r<-0 until n do
        assertEquals(java.lang.Double.doubleToRawLongBits(summary.mean(r)),java.lang.Double.doubleToRawLongBits(expectedMeans(r)))
        for c<-0 until n do
          assertEquals(java.lang.Double.doubleToRawLongBits(summary.covariance(r,c)),java.lang.Double.doubleToRawLongBits(expected(r)(c)))
      val reverse=ok(CenteredMoments.make(n));samples.reverse.foreach((v,w)=>ok(reverse.add(v,w)))
      val a=ok(CenteredMoments.make(n));val b=ok(CenteredMoments.make(n))
      samples.zipWithIndex.foreach: (p,i)=>
        ok((if i%2==0 then a else b).add(p._1,p._2))
      ok(a.merge(b))
      for candidate<-Vector(ok(reverse.finish(0)),ok(a.finish(0))) do
        for r<-0 until n;c<-0 until n do
          assertEquals(java.lang.Double.doubleToRawLongBits(candidate.covariance(r,c)),java.lang.Double.doubleToRawLongBits(summary.covariance(r,c)))

  test("fixed limbs span every finite exponent, both signs and extreme weighted carries"):
    for exponent<-0 to 2046 by 17 do
      val positive=java.lang.Double.longBitsToDouble((exponent.toLong<<52)|0xfffffffffffffL)
      val negative=java.lang.Double.longBitsToDouble((exponent.toLong<<52)|0x8000000000000001L)
      verify(Vector((Array(positive,negative),1L),(Array(negative,positive),3L),
        (Array(java.lang.Double.MIN_VALUE,-java.lang.Double.MIN_VALUE),7L)))
    verify(Vector((Array(Double.MaxValue,-Double.MaxValue),CenteredMoments.MaxFrequency)))
    verify(Vector((Array(1e100,-1e-100),CenteredMoments.MaxFrequency-2),(Array(-1e100,1e-100),2L)))
    verify(Vector((Array(0.0,-0.0),CenteredMoments.MaxFrequency)))

  test("mixed exponents cancellation and sparse carries agree with full-unit oracle"):
    val random=new scala.util.Random(81)
    val samples=Vector.tabulate(100): _=>
      val values=Array.tabulate(3): _=>
        var bits=random.nextLong()
        while ((bits>>>52)&0x7ffL)==0x7ffL do bits=random.nextLong()
        java.lang.Double.longBitsToDouble(bits)
      (values,1L+(random.nextInt(10000).toLong))
    verify(samples)
    verify(Vector.tabulate(40)(i=>(Array(java.lang.Math.scalb(1.0,i-20),-java.lang.Math.scalb(1.0,i-20)),1L<<i)))

  test("chunked physical state crosses cell1024 with independent ownership and exact merge"):
    verify(Vector.tabulate(4)(t=>(Array.tabulate(65)(r=>java.lang.Math.scalb((r+1).toDouble*(t-1),-30)),t+1L)))

  test("failed and zero-frequency adds do not alter moments; copied and merged source state is independent"):
    val state=ok(CenteredMoments.make(2));val values=Array(1.0,-2.0)
    ok(state.add(values,3));val before=ok(state.finish(0));val copy=state.copy()
    assert(state.add(Array(8.0,Double.NaN),1).isLeft)
    assert(state.add(values,CenteredMoments.MaxFrequency).isLeft)
    ok(state.add(Array(9.0,10.0),0));values(0)=99
    assertEquals(state.frequency,3L);assertEquals(state.positiveOccurrences,1L)
    assertEquals(ok(state.finish(0)).mean(0),before.mean(0))
    ok(copy.add(Array(3.0,4.0)));assertEquals(state.frequency,3L)
    val merged=state.copy();ok(merged.merge(copy));assertEquals(copy.frequency,4L)
    assert(merged.merge(merged).isLeft)

  test("resources bound padded signed arrays, retained scratch and copy; original logical shape boundary is preserved"):
    val resource=ok(CenteredMoments.resources(2))
    assert(resource.stateBytes>=BigInt(2*72+3*142)*4+2*13+68+16)
    assertEquals(resource.copyAdditionalBytes,resource.stateBytes)
    assertEquals(resource.addScratchBytes,BigInt(0));assertEquals(resource.mergeScratchBytes,BigInt(0))
    assert(CenteredMoments.resources(46340).isRight)
    assert(CenteredMoments.resources(46341).isLeft)
