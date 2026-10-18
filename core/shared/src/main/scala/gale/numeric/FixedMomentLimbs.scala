package gale.numeric

/** Fixed-width radix2^30 two's-complement cells. Spare sign bits contain every
  * admitted exact sum/product; modular carries give exact signed cancellation.
  * Every finite moment is bounded by admitted total frequency, so carries cannot
  * leave a cell. Chunks keep each physical array inside Int storage for every
  * dimension already admitted by CenteredMoments' dense/triangular shape checks.
  */
private[numeric] final class FixedMomentLimbs private (
    val cells: Int,val digits: Int,private val chunks: Array[Array[Int]]):
  import FixedMomentLimbs.*
  def addShifted(cell: Int,term: Array[Int],length: Int,shift: Int,negative: Boolean): Unit =
    if length>0 then
      val data=chunks(cell/CellsPerChunk)
      val base=(cell%CellsPerChunk)*digits
      val first=shift/Bits
      val offset=shift%Bits
      var i=0;var carry=0L;var shiftedCarry=0L
      while i<length do
        val shifted=term(i).toLong<<offset
        val digit=(shifted&Mask)+shiftedCarry
        val sum=data(base+first+i).toLong+(if negative then -digit else digit)+carry
        data(base+first+i)=(sum&Mask).toInt
        shiftedCarry=shifted>>>Bits;carry=sum>>Bits
        i+=1
      var index=base+first+i
      while index<base+digits && (carry!=0 || shiftedCarry!=0) do
        val sum=data(index).toLong+(if negative then -shiftedCarry else shiftedCarry)+carry
        data(index)=(sum&Mask).toInt;carry=sum>>Bits;shiftedCarry=0;index+=1
      // Discard unsigned carry/borrow at the fixed width. The mathematical
      // frequency bound proves no signed overflow: this is exact modulo2^B.

  def merge(other: FixedMomentLimbs): Unit =
    var cell=0
    while cell<cells do
      val a=chunks(cell/CellsPerChunk);val b=other.chunks(cell/CellsPerChunk)
      val base=(cell%CellsPerChunk)*digits
      var i=0;var carry=0L
      while i<digits do
        val sum=a(base+i).toLong+b(base+i).toLong+carry
        a(base+i)=(sum&Mask).toInt;carry=sum>>>Bits;i+=1
      cell+=1
  def copy(): FixedMomentLimbs = new FixedMomentLimbs(cells,digits,chunks.map(_.clone()))
  def integer(cell: Int): BigInt =
    val data=chunks(cell/CellsPerChunk);val base=(cell%CellsPerChunk)*digits
    var i=digits-1
    while i>=0 && data(base+i)==0 do i-=1
    var result=BigInt(0)
    while i>=0 do
      result=(result<<Bits)+data(base+i);i-=1
    if (data(base+digits-1)&(1<<29))!=0 then result-(BigInt(1)<<(digits*Bits)) else result

private[numeric] object FixedMomentLimbs:
  val Bits=30
  val Mask: Long=(1L<<Bits)-1
  private val CellsPerChunk=1024
  def zero(cells: Int,digits: Int): FixedMomentLimbs =
    val count=((cells.toLong+CellsPerChunk-1)/CellsPerChunk).toInt
    val chunks=Array.tabulate(count)(i=>new Array[Int](math.min(CellsPerChunk.toLong,cells.toLong-i.toLong*CellsPerChunk).toInt*digits))
    new FixedMomentLimbs(cells,digits,chunks)
  def length(values: Array[Int]): Int =
    var n=values.length
    while n>0 && values(n-1)==0 do n-=1
    n
  /** Unsigned radix30 products fit a signed Long: (2^30-1)^2 plus
    * an output digit and incoming carry stays below2^61. JVM and JS agree.
    */
  def multiply(a: Array[Int],aLength: Int,b: Array[Int],bLength: Int,out: Array[Int]): Int =
    java.util.Arrays.fill(out,0)
    var i=0
    while i<aLength do
      var j=0;var carry=0L
      while j<bLength do
        val product=a(i).toLong*b(j).toLong+out(i+j).toLong+carry
        out(i+j)=(product&Mask).toInt;carry=product>>>Bits;j+=1
      var index=i+bLength
      while carry!=0 do
        val sum=out(index).toLong+carry
        out(index)=(sum&Mask).toInt;carry=sum>>>Bits;index+=1
      i+=1
    length(out)
