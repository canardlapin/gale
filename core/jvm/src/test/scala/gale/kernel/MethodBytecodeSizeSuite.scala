package gale.kernel

import java.io.ByteArrayInputStream
import java.io.DataInputStream
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import scala.jdk.CollectionConverters.*

/** Guards every compiled core method against HotSpot's huge-method limit.
  *
  * HotSpot never JIT-compiles a method whose bytecode exceeds 8000 bytes
  * (`-XX:+DontCompileHugeMethods`, the default), so such a method runs in the
  * interpreter for the life of the JVM. A heavily inlined numerical kernel can
  * silently cross the limit and become 10-100x slower; this suite reads the
  * `Code` attribute of every method in the core class directory and fails on any
  * method over the limit.
  */
class MethodBytecodeSizeSuite extends munit.FunSuite:
  private val HugeMethodLimit = 8000

  private def classRoot: Path =
    Paths.get(classOf[gale.linalg.LinAlgError].getProtectionDomain.getCodeSource.getLocation.toURI)

  /** `(methodName, codeLength)` for every method with a `Code` attribute. */
  private def methodCodeLengths(bytes: Array[Byte]): List[(String, Int)] =
    val in = new DataInputStream(new ByteArrayInputStream(bytes))
    require(in.readInt() == 0xcafebabe, "not a class file")
    in.readUnsignedShort()
    in.readUnsignedShort()
    val cpCount = in.readUnsignedShort()
    val utf8 = new Array[String](cpCount)
    var i = 1
    while i < cpCount do
      in.readUnsignedByte() match
        case 1 => utf8(i) = in.readUTF()
        case 3 | 4 => in.skipBytes(4)
        case 5 | 6 =>
          in.skipBytes(8)
          i += 1
        case 7 | 8 | 16 | 19 | 20 => in.skipBytes(2)
        case 9 | 10 | 11 | 12 | 17 | 18 => in.skipBytes(4)
        case 15 => in.skipBytes(3)
        case tag => throw new IllegalStateException(s"unknown constant-pool tag $tag")
      i += 1
    in.skipBytes(6) // access flags, this class, super class
    in.skipBytes(2 * in.readUnsignedShort()) // interfaces
    def skipAttributes(): Unit =
      val count = in.readUnsignedShort()
      var a = 0
      while a < count do
        in.skipBytes(2)
        in.skipBytes(in.readInt())
        a += 1
    val fieldCount = in.readUnsignedShort()
    var f = 0
    while f < fieldCount do
      in.skipBytes(6)
      skipAttributes()
      f += 1
    val methodCount = in.readUnsignedShort()
    val out = List.newBuilder[(String, Int)]
    var m = 0
    while m < methodCount do
      in.skipBytes(2)
      val name = utf8(in.readUnsignedShort())
      in.skipBytes(2)
      val attrCount = in.readUnsignedShort()
      var a = 0
      while a < attrCount do
        val attrName = utf8(in.readUnsignedShort())
        val len = in.readInt()
        if attrName == "Code" then
          in.skipBytes(4) // max stack, max locals
          out += name -> in.readInt()
          in.skipBytes(len - 8)
        else in.skipBytes(len)
        a += 1
      m += 1
    out.result()

  test("no compiled core method exceeds the HotSpot huge-method limit"):
    val root = classRoot
    assume(Files.isDirectory(root), s"core classes are not a directory: $root")
    val classFiles = Files.walk(root).iterator.asScala.filter(_.toString.endsWith(".class")).toList
    assert(classFiles.size > 100, s"expected the core class directory, found ${classFiles.size} classes")
    val huge =
      for
        file <- classFiles
        (method, length) <- methodCodeLengths(Files.readAllBytes(file))
        if length > HugeMethodLimit
      yield s"${root.relativize(file)}#$method: $length bytes"
    assert(huge.isEmpty, s"methods over $HugeMethodLimit bytes are never JIT-compiled:\n${huge.mkString("\n")}")

  test("the class-file reader measures a known method"):
    val bytes = Files.readAllBytes(classRoot.resolve("gale/spectral/DenseSpectralKernels$.class"))
    val lengths = methodCodeLengths(bytes).toMap
    assert(lengths.get("doubleQrSweep").exists(_ > 500), lengths.keys.toString)
