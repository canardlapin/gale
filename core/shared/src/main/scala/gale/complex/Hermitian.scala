package gale.complex

/** Explicit cutoff in the input matrix's coordinates, with no ridge or hidden floor. */
enum HermitianCutoff:
  case Relative(value: Double)
  case Absolute(value: Double)
final case class HermitianPolicy(cutoff: HermitianCutoff, maxSweeps: Int = 100,
    convergenceTolerance: Double = 1e-13, symmetryTolerance: Double = 1e-12,
    psdTolerance: Double = 1e-12, certificateTolerance: Double = 1e-9)
final case class HermitianRank(dimension: Int, eigenvalues: Vector[Double], retained: Vector[Int],
    cutoff: HermitianCutoff, threshold: Double, acceptedNegativeRoundoff: Vector[Double]):
  def rank: Int = retained.size
final case class HermitianCertificate(sweeps: Int, offDiagonal: Double, eigenResidual: Double,
    orthogonalityError: Double, relativeHermitianCorrection: Double)
enum HermitianSolveMode:
  case FullRank, RetainedLeastSquares
final class HermitianSolve private[complex] (val values: OwnedComplex, val mode: HermitianSolveMode,
    val relativeResidual: Double, val relativeDiscardedRhs: Double, val rank: HermitianRank, val numericalBytes: BigInt)

/** Certified portable Hermitian eigensystem and PSD spectral factor. */
final class HermitianFactor private[complex] (val input: OwnedComplex, val eigenvectors: OwnedComplex,
    val record: HermitianRank, val certificate: HermitianCertificate, val policy: HermitianPolicy,
    val numericalBytes: BigInt):
  val backend: String = "portable-complex-jacobi-v1"
  private def spectralMatrix(weights: Int => Double, budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    val n = input.rows
    OwnedComplex.inspect(n,n,budget).flatMap: _ =>
      val re = new Array[Double](n*n); val im = new Array[Double](n*n)
      var i = 0; var invalid = false
      while i < n do
        var j = 0
        while j <= i do
          var rr = 0.0; var ii = 0.0
          record.retained.foreach: k =>
            val a = eigenvectors(i,k); val b = eigenvectors(j,k); val w = weights(k)
            rr += (a.re*b.re+a.im*b.im)*w; ii += (a.im*b.re-a.re*b.im)*w
          re(i*n+j) = rr; im(i*n+j) = if i==j then 0.0 else ii
          re(j*n+i) = rr; im(j*n+i) = if i==j then 0.0 else -ii
          invalid ||= !rr.isFinite || !ii.isFinite
          j += 1
        i += 1
      if invalid then Left(ComplexFailure.Numerical("nonfinite spectral matrix"))
      else Right(new OwnedComplex(n,n,re,im))
  def projector(budget: BigInt): Either[ComplexFailure, OwnedComplex] = spectralMatrix(_=>1.0,budget)
  def pseudoInverse(budget: BigInt): Either[ComplexFailure, OwnedComplex] = spectralMatrix(k=>1.0/record.eigenvalues(k),budget)
  /** F F^H equals the retained PSD part; discarded eigenvalues remain in the rank record. */
  def psdFactor(budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    val n = input.rows; val rank = record.rank
    OwnedComplex.inspect(n,rank,budget).map: _ =>
      val re = new Array[Double](n*rank); val im = new Array[Double](n*rank)
      var i=0
      while i<n do
        var j=0
        while j<rank do
          val k=record.retained(j); val scale=math.sqrt(record.eigenvalues(k)); val v=eigenvectors(i,k)
          re(i*rank+j)=v.re*scale; im(i*rank+j)=v.im*scale; j+=1
        i+=1
      new OwnedComplex(n,rank,re,im)
  def solve(rhs: OwnedComplex,mode: HermitianSolveMode,budget: BigInt): Either[ComplexFailure,HermitianSolve] =
    val n=input.rows; val m=rhs.cols; val bytes=BigInt(16)*(4*BigInt(n)*m+BigInt(n)*n)
    if rhs.rows!=n then return Left(ComplexFailure.Invalid("Hermitian RHS rows differ"))
    if bytes>budget then return Left(ComplexFailure.Budget(bytes,budget))
    if mode==HermitianSolveMode.FullRank && record.rank!=n then return Left(ComplexFailure.RankDeficient(record.rank,n))
    // Explicit spectral projection/solve; no regularization of discarded directions.
    val rr=new Array[Double](n*m); val ii=new Array[Double](n*m)
    val pr=new Array[Double](n*m); val pi=new Array[Double](n*m)
    var c=0
    while c<m do
      record.retained.foreach: k =>
        var ar=0.0; var ai=0.0; var i=0
        while i<n do
          val v=eigenvectors(i,k); val b=rhs(i,c)
          ar+=v.re*b.re+v.im*b.im; ai+=v.re*b.im-v.im*b.re; i+=1
        i=0
        while i<n do
          val v=eigenvectors(i,k); val real=v.re*ar-v.im*ai; val imag=v.re*ai+v.im*ar
          pr(i*m+c)+=real; pi(i*m+c)+=imag
          rr(i*m+c)+=real/record.eigenvalues(k); ii(i*m+c)+=imag/record.eigenvalues(k); i+=1
      c+=1
    if !rr.forall(_.isFinite) || !ii.forall(_.isFinite) || !pr.forall(_.isFinite) || !pi.forall(_.isFinite) then
      return Left(ComplexFailure.Numerical("nonfinite Hermitian solve"))
    val x=new OwnedComplex(n,m,rr,ii)
    input.multiply(x,BigInt(n)*m*16).flatMap: actual =>
      var residual=0.0; var discarded=0.0; var certificate=0.0; var scale=0.0; var k=0
      while k<rr.length do
        scale=math.max(scale,math.hypot(rhs.re(k),rhs.im(k)))
        residual=math.max(residual,math.hypot(actual.re(k)-rhs.re(k),actual.im(k)-rhs.im(k)))
        discarded=math.max(discarded,math.hypot(pr(k)-rhs.re(k),pi(k)-rhs.im(k)))
        certificate=math.max(certificate,math.hypot(actual.re(k)-pr(k),actual.im(k)-pi(k))); k+=1
      val denominator=if scale==0 then 1.0 else scale
      if !certificate.isFinite || certificate/denominator>policy.certificateTolerance then
        Left(ComplexFailure.Numerical("Hermitian solve residual certificate failed"))
      else Right(new HermitianSolve(x,mode,residual/denominator,discarded/denominator,record,bytes))

object Hermitian:
  val MaxDimension: Int=128
  /** Conservative primitive peak including projected input, rotations, sorting,
    * certificates and retained owned outputs; excludes caller input and GC/objects.
    */
  def numericalBytes(n: Int): BigInt=BigInt(16)*(16*BigInt(n)*n+8*BigInt(n))
  def inspect(n: Int,p: HermitianPolicy,budget: BigInt): Either[ComplexFailure,BigInt] =
    val cutoff=p.cutoff match
      case HermitianCutoff.Relative(v)=>v
      case HermitianCutoff.Absolute(v)=>v
    if n<1 || n>MaxDimension || p.maxSweeps<1 || !cutoff.isFinite || cutoff<0 ||
      !p.convergenceTolerance.isFinite || p.convergenceTolerance<1e-15 || p.convergenceTolerance>1e-6 ||
      !p.symmetryTolerance.isFinite || p.symmetryTolerance<0 || p.symmetryTolerance>1e-6 ||
      !p.psdTolerance.isFinite || p.psdTolerance<0 || p.psdTolerance>1e-6 ||
      !p.certificateTolerance.isFinite || p.certificateTolerance<=0 || p.certificateTolerance>1e-6 then
      Left(ComplexFailure.Invalid("bounded Hermitian dimension and explicit finite tolerances/cutoff required"))
    else if numericalBytes(n)>budget then Left(ComplexFailure.Budget(numericalBytes(n),budget))
    else Right(numericalBytes(n))
  def factor(matrix: OwnedComplex,policy: HermitianPolicy,budget: BigInt,
      cancelled: ()=>Boolean=()=>false): Either[ComplexFailure,HermitianFactor] =
    val n=matrix.rows
    val bytes=inspect(n,policy,budget) match
      case Left(e)=>return Left(e)
      case Right(v)=>v
    if matrix.cols!=n then return Left(ComplexFailure.Invalid("square Hermitian matrix required"))
    if cancelled() then return Left(ComplexFailure.Cancelled)
    var scale=0.0; var k=0
    while k<n*n do
      scale=math.max(scale,math.hypot(matrix.re(k),matrix.im(k))); k+=1
    if !scale.isFinite then return Left(ComplexFailure.Numerical("Hermitian magnitude overflow"))
    val divisor=if scale==0 then 1.0 else scale
    val hr=new Array[Double](n*n); val hi=new Array[Double](n*n)
    var correction=0.0; var i=0
    while i<n do
      var j=0
      while j<=i do
        val ar=matrix.re(i*n+j)/divisor; val ai=matrix.im(i*n+j)/divisor
        val br=matrix.re(j*n+i)/divisor; val bi=matrix.im(j*n+i)/divisor
        val error=math.hypot(ar-br,ai+bi)
        if error>policy.symmetryTolerance then return Left(ComplexFailure.Invalid("matrix is not Hermitian within declared tolerance"))
        correction=math.max(correction,error/2)
        hr(i*n+j)=(ar+br)/2; hr(j*n+i)=hr(i*n+j)
        hi(i*n+j)=if i==j then 0.0 else (ai-bi)/2; hi(j*n+i) = -hi(i*n+j); j+=1
      i+=1
    val originalR=hr.clone(); val originalI=hi.clone()
    val vr=new Array[Double](n*n); val vi=new Array[Double](n*n)
    i=0; while i<n do {vr(i*n+i)=1.0;i+=1}
    def off: Double =
      var largest=0.0;var a=0
      while a<n do
        var b=a+1
        while b<n do {largest=math.max(largest,math.hypot(hr(a*n+b),hi(a*n+b)));b+=1}
        a+=1
      largest
    var sweeps=0; var largest=off
    while largest>policy.convergenceTolerance && sweeps<policy.maxSweeps do
      var p=0
      while p<n do
        var q=p+1
        while q<n do
          if cancelled() then return Left(ComplexFailure.Cancelled)
          val rho=math.hypot(hr(p*n+q),hi(p*n+q))
          if rho>policy.convergenceTolerance then
            val er=hr(p*n+q)/rho;val ei=hi(p*n+q)/rho
            val tau=(hr(q*n+q)-hr(p*n+p))/(2*rho)
            val t=(if tau>=0 then 1.0 else -1.0)/(math.abs(tau)+math.hypot(1.0,tau))
            val c=1/math.hypot(1.0,t);val s=t*c
            k=0
            while k<n do
              if k!=p && k!=q then
                val ar=hr(k*n+p);val ai=hi(k*n+p);val br=hr(k*n+q);val bi=hi(k*n+q)
                val pr=c*ar-s*(er*br+ei*bi);val pi=c*ai-s*(er*bi-ei*br)
                val qr=s*(er*ar-ei*ai)+c*br;val qi=s*(er*ai+ei*ar)+c*bi
                hr(k*n+p)=pr;hi(k*n+p)=pi;hr(p*n+k)=pr;hi(p*n+k) = -pi
                hr(k*n+q)=qr;hi(k*n+q)=qi;hr(q*n+k)=qr;hi(q*n+k) = -qi
              val ar=vr(k*n+p);val ai=vi(k*n+p);val br=vr(k*n+q);val bi=vi(k*n+q)
              vr(k*n+p)=c*ar-s*(er*br+ei*bi);vi(k*n+p)=c*ai-s*(er*bi-ei*br)
              vr(k*n+q)=s*(er*ar-ei*ai)+c*br;vi(k*n+q)=s*(er*ai+ei*ar)+c*bi
              k+=1
            hr(p*n+p)-=t*rho;hr(q*n+q)+=t*rho
            hr(p*n+q)=0;hi(p*n+q)=0;hr(q*n+p)=0;hi(q*n+p)=0
          q+=1
        p+=1
      sweeps+=1;largest=off
    if largest>policy.convergenceTolerance then return Left(ComplexFailure.NonConverged(sweeps,largest))
    val order=(0 until n).sortBy(k=>hr(k*n+k)).toVector
    val normalized=order.map(k=>hr(k*n+k));val values=normalized.map(_*divisor)
    if !values.forall(_.isFinite) then return Left(ComplexFailure.Numerical("eigenvalue rescaling overflow"))
    val sr=new Array[Double](n*n);val si=new Array[Double](n*n)
    i=0
    while i<n do
      k=0
      while k<n do {sr(i*n+k)=vr(i*n+order(k));si(i*n+k)=vi(i*n+order(k));k+=1}
      i+=1
    var eigenError=0.0;var orthogonality=0.0
    i=0
    while i<n do
      if cancelled() then return Left(ComplexFailure.Cancelled)
      var j=0
      while j<n do
        var ar=0.0;var ai=0.0;var gr=0.0;var gi=0.0;k=0
        while k<n do
          ar+=originalR(i*n+k)*sr(k*n+j)-originalI(i*n+k)*si(k*n+j)
          ai+=originalR(i*n+k)*si(k*n+j)+originalI(i*n+k)*sr(k*n+j)
          gr+=sr(k*n+i)*sr(k*n+j)+si(k*n+i)*si(k*n+j)
          gi+=sr(k*n+i)*si(k*n+j)-si(k*n+i)*sr(k*n+j);k+=1
        eigenError=math.max(eigenError,math.hypot(ar-sr(i*n+j)*normalized(j),ai-si(i*n+j)*normalized(j)))
        orthogonality=math.max(orthogonality,math.hypot(gr-(if i==j then 1.0 else 0.0),gi));j+=1
      i+=1
    if !eigenError.isFinite || !orthogonality.isFinite || eigenError>policy.certificateTolerance || orthogonality>policy.certificateTolerance then
      return Left(ComplexFailure.Numerical("Hermitian eigensystem residual or orthogonality certificate failed"))
    val max=values.last.max(0.0);val negative=policy.psdTolerance*max
    if values.head < -negative then return Left(ComplexFailure.Indefinite(values.head,negative))
    val threshold=policy.cutoff match
      case HermitianCutoff.Relative(v)=>v*max
      case HermitianCutoff.Absolute(v)=>v
    if !threshold.isFinite then return Left(ComplexFailure.Numerical("rank threshold overflow"))
    val kept=values.indices.filter(k=>values(k)>threshold).toVector
    val record=HermitianRank(n,values,kept,policy.cutoff,threshold,values.filter(_<0))
    val cert=HermitianCertificate(sweeps,largest,eigenError,orthogonality,correction)
    if cancelled() then return Left(ComplexFailure.Cancelled)
    // Keep the supplied input, not the Hermitian projection, for solve residual checks.
    Right(new HermitianFactor(matrix,new OwnedComplex(n,n,sr,si),record,cert,policy,bytes))
