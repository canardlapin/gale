package gale.backend.jvm.blas;

import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.foreign.ValueLayout;
import java.lang.invoke.MethodHandle;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

/** Minimal, JDK-22-final FFM bindings for the three coarse CBLAS seams. */
public final class CblasBindings implements AutoCloseable {
  public static final int ROW_MAJOR = 101;
  public static final int COL_MAJOR = 102;
  public static final int NO_TRANS = 111;
  public static final int TRANS = 112;
  public static final int UPPER = 121;

  private static final FunctionDescriptor DGEMM = FunctionDescriptor.ofVoid(
      ValueLayout.JAVA_INT, ValueLayout.JAVA_INT, ValueLayout.JAVA_INT,
      ValueLayout.JAVA_INT, ValueLayout.JAVA_INT, ValueLayout.JAVA_INT,
      ValueLayout.JAVA_DOUBLE, ValueLayout.ADDRESS, ValueLayout.JAVA_INT,
      ValueLayout.ADDRESS, ValueLayout.JAVA_INT, ValueLayout.JAVA_DOUBLE,
      ValueLayout.ADDRESS, ValueLayout.JAVA_INT);
  private static final FunctionDescriptor DGEMV = FunctionDescriptor.ofVoid(
      ValueLayout.JAVA_INT, ValueLayout.JAVA_INT, ValueLayout.JAVA_INT,
      ValueLayout.JAVA_INT, ValueLayout.JAVA_DOUBLE, ValueLayout.ADDRESS,
      ValueLayout.JAVA_INT, ValueLayout.ADDRESS, ValueLayout.JAVA_INT,
      ValueLayout.JAVA_DOUBLE, ValueLayout.ADDRESS, ValueLayout.JAVA_INT);
  private static final FunctionDescriptor DSYRK = FunctionDescriptor.ofVoid(
      ValueLayout.JAVA_INT, ValueLayout.JAVA_INT, ValueLayout.JAVA_INT,
      ValueLayout.JAVA_INT, ValueLayout.JAVA_INT, ValueLayout.JAVA_DOUBLE,
      ValueLayout.ADDRESS, ValueLayout.JAVA_INT, ValueLayout.JAVA_DOUBLE,
      ValueLayout.ADDRESS, ValueLayout.JAVA_INT);
  private static final FunctionDescriptor SET_THREADS =
      FunctionDescriptor.ofVoid(ValueLayout.JAVA_INT);

  private final String libraryName;
  private final Arena libraryArena;
  private final MethodHandle dgemm;
  private final MethodHandle dgemv;
  private final MethodHandle dsyrk;
  private final MethodHandle setThreads;
  private final String threadSetterName;
  private static Integer configuredThreadCount;

  private CblasBindings(String libraryName, Arena arena, SymbolLookup lookup) {
    Linker linker = Linker.nativeLinker();
    this.libraryName = libraryName;
    this.libraryArena = arena;
    this.dgemm = linker.downcallHandle(required(lookup, "cblas_dgemm"), DGEMM);
    this.dgemv = linker.downcallHandle(required(lookup, "cblas_dgemv"), DGEMV);
    this.dsyrk = linker.downcallHandle(required(lookup, "cblas_dsyrk"), DSYRK);
    String[] setters = {"openblas_set_num_threads", "goto_set_num_threads", "MKL_Set_Num_Threads", "mkl_set_num_threads"};
    MethodHandle found = null;
    String foundName = null;
    for (String name : setters) {
      Optional<MemorySegment> symbol = lookup.find(name);
      if (symbol.isPresent()) {
        found = linker.downcallHandle(symbol.get(), SET_THREADS);
        foundName = name;
        break;
      }
    }
    this.setThreads = found;
    this.threadSetterName = foundName;
  }

  public static CblasBindings loadDefault() {
    return load(candidates());
  }

  public static CblasBindings load(List<String> requestedCandidates) {
    List<String> failures = new ArrayList<>();
    for (String candidate : new LinkedHashSet<>(requestedCandidates)) {
      if (candidate == null || candidate.isBlank()) continue;
      Arena arena = Arena.ofShared();
      try {
        SymbolLookup lookup = SymbolLookup.libraryLookup(candidate, arena);
        return new CblasBindings(candidate, arena, lookup);
      } catch (Throwable error) {
        failures.add(candidate + ": " + error.getClass().getSimpleName() + ": " + error.getMessage());
        if (arena.scope().isAlive()) arena.close();
      }
    }
    throw new IllegalStateException(
        "No usable CBLAS library found. Set -Dgale.blas.library=/absolute/path/to/library. " +
        "Probe failures: " + String.join(" | ", failures));
  }

  public String libraryName() { return libraryName; }
  public boolean hasThreadControl() { return setThreads != null; }
  public String threadSetterName() { return threadSetterName; }

  public boolean configureThreads(int threads) {
    if (setThreads == null) return false;
    synchronized (CblasBindings.class) {
      if (configuredThreadCount != null && configuredThreadCount != threads) {
        throw new IllegalStateException(
            "native BLAS thread count is process-global and was already configured to " +
            configuredThreadCount + ", cannot change it to " + threads);
      }
    }
    try {
      setThreads.invokeExact(threads);
      synchronized (CblasBindings.class) {
        configuredThreadCount = threads;
      }
      return true;
    } catch (Throwable error) {
      throw failure("set native thread count", error);
    }
  }

  public void dgemm(int layout, int transA, int transB, int m, int n, int k, double alpha,
                    MemorySegment a, int lda, MemorySegment b, int ldb,
                    double beta, MemorySegment c, int ldc) {
    try {
      dgemm.invokeExact(layout, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } catch (Throwable error) {
      throw failure("cblas_dgemm", error);
    }
  }

  public void dgemv(int trans, int m, int n, double alpha, MemorySegment a, int lda,
                    MemorySegment x, int incX, double beta, MemorySegment y, int incY) {
    try {
      dgemv.invokeExact(ROW_MAJOR, trans, m, n, alpha, a, lda, x, incX, beta, y, incY);
    } catch (Throwable error) {
      throw failure("cblas_dgemv", error);
    }
  }

  public void dsyrk(int n, int k, double alpha, MemorySegment a, int lda,
                    double beta, MemorySegment c, int ldc) {
    try {
      dsyrk.invokeExact(ROW_MAJOR, UPPER, TRANS, n, k, alpha, a, lda, beta, c, ldc);
    } catch (Throwable error) {
      throw failure("cblas_dsyrk", error);
    }
  }

  @Override public void close() {
    if (libraryArena.scope().isAlive()) libraryArena.close();
  }

  private static MemorySegment required(SymbolLookup lookup, String name) {
    return lookup.find(name).orElseThrow(() -> new IllegalArgumentException("missing required symbol " + name));
  }

  private static RuntimeException failure(String operation, Throwable cause) {
    return new IllegalStateException(operation + " failed through FFM", cause);
  }

  private static List<String> candidates() {
    LinkedHashSet<String> out = new LinkedHashSet<>();
    add(out, System.getProperty("gale.blas.library"));
    add(out, System.getenv("GALE_BLAS_LIBRARY"));
    String configured = System.getProperty("gale.blas.candidates");
    if (configured != null) for (String item : configured.split(",")) add(out, item);

    String os = System.getProperty("os.name", "").toLowerCase(Locale.ROOT);
    if (os.contains("mac")) {
      add(out, "/opt/homebrew/opt/openblas/lib/libopenblas.dylib");
      add(out, "/usr/local/opt/openblas/lib/libopenblas.dylib");
      add(out, "/System/Library/Frameworks/Accelerate.framework/Accelerate");
      add(out, "/usr/lib/libblas.dylib");
      add(out, "libopenblas.dylib");
      add(out, "libmkl_rt.dylib");
    } else if (os.contains("win")) {
      add(out, "libopenblas.dll");
      add(out, "openblas.dll");
      add(out, "mkl_rt.dll");
    } else {
      add(out, "/usr/lib/x86_64-linux-gnu/libopenblas.so.0");
      add(out, "/usr/lib/aarch64-linux-gnu/libopenblas.so.0");
      add(out, "/usr/lib64/libopenblas.so.0");
      add(out, "/usr/lib/libopenblas.so.0");
      add(out, "libopenblas.so.0");
      add(out, "libblas.so.3");
      add(out, "libmkl_rt.so");
    }
    return List.copyOf(out);
  }

  private static void add(LinkedHashSet<String> out, String value) {
    if (value != null && !value.isBlank()) out.add(value.trim());
  }
}
