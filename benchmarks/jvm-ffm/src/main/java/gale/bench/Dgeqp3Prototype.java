package gale.bench;

import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.foreign.ValueLayout;
import java.lang.invoke.MethodHandle;

/** Benchmark-only FFM binding for LAPACK dgeqp3. It is deliberately isolated
 * from Gale's production backend and owns no memory beyond one factor call.
 */
public final class Dgeqp3Prototype implements AutoCloseable {
  private static final String ACCELERATE =
      "/System/Library/Frameworks/Accelerate.framework/Accelerate";
  private static final FunctionDescriptor DGEQP3 = FunctionDescriptor.ofVoid(
      ValueLayout.ADDRESS, ValueLayout.ADDRESS, ValueLayout.ADDRESS,
      ValueLayout.ADDRESS, ValueLayout.ADDRESS, ValueLayout.ADDRESS,
      ValueLayout.ADDRESS, ValueLayout.ADDRESS, ValueLayout.ADDRESS);

  public record Result(
      double[] packedColumnMajor,
      double[] tau,
      int[] permutation,
      int lwork,
      long peakNativeBytes) {}

  private final String libraryName;
  private final Arena libraryArena;
  private final MethodHandle dgeqp3;

  private Dgeqp3Prototype(String libraryName, Arena libraryArena, MethodHandle dgeqp3) {
    this.libraryName = libraryName;
    this.libraryArena = libraryArena;
    this.dgeqp3 = dgeqp3;
  }

  public static Dgeqp3Prototype loadDefault() {
    String configured = System.getProperty("gale.blas.library");
    String library = configured == null || configured.isBlank() ? ACCELERATE : configured.trim();
    Arena arena = Arena.ofShared();
    try {
      SymbolLookup lookup = SymbolLookup.libraryLookup(library, arena);
      MemorySegment symbol = lookup.find("dgeqp3_")
          .or(() -> lookup.find("dgeqp3"))
          .or(() -> lookup.find("DGEQP3"))
          .orElseThrow(() -> new IllegalArgumentException("missing LAPACK dgeqp3 symbol in " + library));
      MethodHandle handle = Linker.nativeLinker().downcallHandle(symbol, DGEQP3);
      return new Dgeqp3Prototype(library, arena, handle);
    } catch (Throwable error) {
      if (arena.scope().isAlive()) arena.close();
      throw new IllegalStateException("could not load dgeqp3 from " + library, error);
    }
  }

  public String libraryName() {
    return libraryName;
  }

  /** Takes a caller-owned column-major copy and includes heap/native copies,
   * LAPACK workspace query and factorization, and owned heap copy-out arrays.
   */
  public Result factorOwnedColumnMajor(double[] packed, int rows, int cols) {
    if (rows < 0 || cols < 0) throw new IllegalArgumentException("negative matrix extent");
    if ((long) rows * cols != packed.length) throw new IllegalArgumentException("matrix storage mismatch");
    int limit = Math.min(rows, cols);
    double[] tau = new double[limit];
    int[] permutation = new int[cols];
    if (rows == 0 || cols == 0) {
      for (int j = 0; j < cols; j++) permutation[j] = j;
      return new Result(packed, tau, permutation, 1, 0L);
    }

    try (Arena arena = Arena.ofConfined()) {
      MemorySegment matrix = doubles(arena, packed.length);
      matrix.copyFrom(MemorySegment.ofArray(packed));
      MemorySegment jpvt = ints(arena, cols);
      MemorySegment tauSegment = doubles(arena, limit);
      MemorySegment query = doubles(arena, 1);

      int queryInfo = invoke(rows, cols, matrix, rows, jpvt, tauSegment, query, -1);
      if (queryInfo != 0) throw new IllegalArgumentException("dgeqp3 workspace query info=" + queryInfo);
      double queriedWork = query.get(ValueLayout.JAVA_DOUBLE, 0L);
      if (!Double.isFinite(queriedWork) || queriedWork < 1.0 || queriedWork > Integer.MAX_VALUE) {
        throw new IllegalArgumentException("invalid dgeqp3 workspace size " + queriedWork);
      }
      int lwork = (int) Math.ceil(queriedWork);
      MemorySegment work = doubles(arena, lwork);

      // Query calls are not specified to preserve A or JPVT. Restore both so
      // the timed operation has one deterministic factorization input.
      matrix.copyFrom(MemorySegment.ofArray(packed));
      jpvt.fill((byte) 0);
      int info = invoke(rows, cols, matrix, rows, jpvt, tauSegment, work, lwork);
      if (info < 0) throw new IllegalArgumentException("dgeqp3 rejected argument " + (-info));
      if (info > 0) throw new IllegalStateException("dgeqp3 failed with info=" + info);

      MemorySegment.ofArray(packed).copyFrom(matrix);
      MemorySegment.ofArray(tau).copyFrom(tauSegment);
      for (int j = 0; j < cols; j++) {
        permutation[j] = jpvt.getAtIndex(ValueLayout.JAVA_INT, j) - 1;
      }
      long peakNativeBytes =
          8L * packed.length + 4L * cols + 8L * limit + 8L + 8L * lwork + 20L;
      return new Result(packed, tau, permutation, lwork, peakNativeBytes);
    }
  }

  private int invoke(
      int rows,
      int cols,
      MemorySegment matrix,
      int leadingDimension,
      MemorySegment jpvt,
      MemorySegment tau,
      MemorySegment work,
      int lwork) {
    try (Arena call = Arena.ofConfined()) {
      MemorySegment m = integer(call, rows);
      MemorySegment n = integer(call, cols);
      MemorySegment lda = integer(call, leadingDimension);
      MemorySegment workSize = integer(call, lwork);
      MemorySegment info = integer(call, 0);
      dgeqp3.invokeExact(m, n, matrix, lda, jpvt, tau, work, workSize, info);
      return info.get(ValueLayout.JAVA_INT, 0L);
    } catch (Throwable error) {
      throw new IllegalStateException("dgeqp3 failed through FFM", error);
    }
  }

  private static MemorySegment doubles(Arena arena, int length) {
    return arena.allocate(Math.max(1L, 8L * length), 8L);
  }

  private static MemorySegment ints(Arena arena, int length) {
    return arena.allocate(Math.max(1L, 4L * length), 4L);
  }

  private static MemorySegment integer(Arena arena, int value) {
    MemorySegment result = arena.allocate(ValueLayout.JAVA_INT);
    result.set(ValueLayout.JAVA_INT, 0L, value);
    return result;
  }

  @Override
  public void close() {
    if (libraryArena.scope().isAlive()) libraryArena.close();
  }
}
