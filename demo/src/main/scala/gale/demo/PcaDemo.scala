package gale.demo

import org.scalajs.dom
import org.scalajs.dom.{document, html, window}

import gale.linalg.{DMat, DVec, Matrix}
import gale.spectral.{Eigen, EigenSelection, EigenVectors}

/** A live, in-browser PCA demo built on gale's public API.
  *
  * Everything numeric here is real gale computation running as Scala.js: the
  * synthetic 5-D dataset is generated in the page, centered with
  * `Matrix.tabulate`, the scatter matrix is formed via `centered.t * centered`,
  * the spectrum comes from `Eigen.eigSymmetric`, and the top-2 projection is a
  * gale matrix-matrix product. The "Re-sample" button regenerates the data and
  * re-runs the whole pipeline, timing it with `performance.now()`.
  *
  * Data generation is deterministic (a hand-rolled seeded LCG + Box-Muller, no
  * `scala.util.Random`), mirroring the repo's determinism conventions.
  */
object PcaDemo:

  // ---------------------------------------------------------------------------
  // Deterministic random numbers: a 64-bit LCG (Knuth's MMIX multiplier) with
  // Box-Muller on top for Gaussians. Seeded, portable, dependency-free.
  // ---------------------------------------------------------------------------

  private final class Lcg(seed: Long):
    private var state: Long = seed * 0x9e3779b97f4a7c15L + 0x2545f4914f6cdd1dL
    private var spare: Double = Double.NaN

    /** Uniform in [0, 1): the top 53 bits of the LCG state. */
    def nextDouble(): Double =
      state = state * 6364136223846793005L + 1442695040888963407L
      (state >>> 11).toDouble / 9007199254740992.0 // 2^53

    /** Standard Gaussian via Box-Muller (with the spare cached). */
    def nextGaussian(): Double =
      if !spare.isNaN then
        val out = spare
        spare = Double.NaN
        out
      else
        val u1 = math.max(nextDouble(), 1e-12)
        val u2 = nextDouble()
        val r = math.sqrt(-2.0 * math.log(u1))
        spare = r * math.sin(2.0 * math.Pi * u2)
        r * math.cos(2.0 * math.Pi * u2)

  // ---------------------------------------------------------------------------
  // Synthetic dataset: three Gaussian clusters in 5-D.
  // ---------------------------------------------------------------------------

  private val dims = 5
  private val perCluster = 60
  private val clusterNames = Array("Cluster A", "Cluster B", "Cluster C")

  private val clusterMeans: Array[Array[Double]] = Array(
    Array(3.0, 0.5, 1.0, -1.5, 0.5),
    Array(-2.5, 3.0, -1.0, 1.0, -0.5),
    Array(0.0, -3.0, 2.0, 2.5, -1.5)
  )
  private val clusterSigmas: Array[Array[Double]] = Array(
    Array(0.9, 0.7, 1.1, 0.6, 0.8),
    Array(0.7, 1.0, 0.6, 0.9, 0.7),
    Array(1.0, 0.8, 0.7, 1.1, 0.9)
  )

  private final case class Dataset(points: DMat, labels: Array[Int]):
    def n: Int = labels.length

  private def generate(seed: Long): Dataset =
    val rng = new Lcg(seed)
    val k = clusterMeans.length
    val n = k * perCluster
    val raw = Array.ofDim[Double](n, dims)
    val labels = new Array[Int](n)
    var row = 0
    var c = 0
    while c < k do
      var i = 0
      while i < perCluster do
        var d = 0
        while d < dims do
          raw(row)(d) = clusterMeans(c)(d) + clusterSigmas(c)(d) * rng.nextGaussian()
          d += 1
        labels(row) = c
        row += 1
        i += 1
      c += 1
    Dataset(Matrix.tabulate(n, dims)((i, j) => raw(i)(j)), labels)

  // ---------------------------------------------------------------------------
  // PCA via gale: center -> scatter (Xt X) -> symmetric eigendecomposition ->
  // project onto the top-2 eigenvectors. Same pipeline as the worked PCA
  // example in docs/examples.md / WorkedExamplesSuite.
  // ---------------------------------------------------------------------------

  private final case class PcaResult(
      projected: DMat, // n x 2: PC1 in column 0, PC2 in column 1
      varExplained1: Double,
      varExplained2: Double,
      elapsedMs: Double
  )

  private def runPca(data: Dataset): Either[String, PcaResult] =
    val t0 = window.performance.now()
    val x = data.points
    val n = x.rows

    val colMeans = Array.tabulate(dims)(j => (0 until n).map(x(_, j)).sum / n)
    val centered = Matrix.tabulate(n, dims)((i, j) => x(i, j) - colMeans(j))

    // Unnormalized scatter matrix: a positive multiple of the sample covariance,
    // hence identical eigenvectors and variance ratios.
    val scatter = centered.t * centered

    Eigen.eigSymmetric(scatter, EigenSelection.All, EigenVectors.Right) match
      case Left(err) => Left(err.toString)
      case Right(eig) =>
        // Eigenvalues come back ascending: the top two live at the last indices.
        val last = eig.eigenvalues.length - 1
        val total = (0 to last).map(i => math.max(eig.eigenvalues(i), 0.0)).sum
        val basis = Matrix.tabulate(dims, 2): (i, j) =>
          eig.eigenvectors(i, last - j)
        val projected = centered * basis
        val elapsed = window.performance.now() - t0
        Right(
          PcaResult(
            projected,
            if total > 0 then math.max(eig.eigenvalues(last), 0.0) / total else 0.0,
            if total > 0 then math.max(eig.eigenvalues(last - 1), 0.0) / total else 0.0,
            elapsed
          )
        )

  // ---------------------------------------------------------------------------
  // Theme (light/dark) — colors from a CVD-validated categorical palette; the
  // page CSS mirrors the same surfaces via prefers-color-scheme.
  // ---------------------------------------------------------------------------

  private final case class Theme(
      surface: String,
      textPrimary: String,
      textSecondary: String,
      grid: String,
      series: Array[String]
  )

  private val lightTheme = Theme(
    surface = "#fcfcfb",
    textPrimary = "#0b0b0b",
    textSecondary = "#52514e",
    grid = "#e4e3df",
    series = Array("#2a78d6", "#008300", "#e87ba4")
  )
  private val darkTheme = Theme(
    surface = "#1a1a19",
    textPrimary = "#ffffff",
    textSecondary = "#c3c2b7",
    grid = "#333331",
    series = Array("#3987e5", "#008300", "#d55181")
  )

  private val darkQuery = window.matchMedia("(prefers-color-scheme: dark)")
  private def theme: Theme = if darkQuery.matches then darkTheme else lightTheme

  // ---------------------------------------------------------------------------
  // Rendering: scatter of PC1/PC2 on a canvas with recessive grid + axes,
  // per-cluster colors, and direct cluster labels at the centroids.
  // ---------------------------------------------------------------------------

  private val canvasWidth = 720
  private val canvasHeight = 460

  private def niceStep(raw: Double): Double =
    val mag = math.pow(10.0, math.floor(math.log10(raw)))
    val norm = raw / mag
    val base = if norm < 1.5 then 1.0 else if norm < 3.5 then 2.0 else if norm < 7.5 then 5.0 else 10.0
    base * mag

  private def render(canvas: html.Canvas, data: Dataset, result: PcaResult): Unit =
    val th = theme
    val dpr = math.max(window.devicePixelRatio, 1.0)
    canvas.width = (canvasWidth * dpr).toInt
    canvas.height = (canvasHeight * dpr).toInt
    val ctx = canvas.getContext("2d").asInstanceOf[dom.CanvasRenderingContext2D]
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)

    ctx.fillStyle = th.surface
    ctx.fillRect(0, 0, canvasWidth.toDouble, canvasHeight.toDouble)

    val n = data.n
    val xs = Array.tabulate(n)(i => result.projected(i, 0))
    val ys = Array.tabulate(n)(i => result.projected(i, 1))
    val padFrac = 0.08
    def bounds(v: Array[Double]): (Double, Double) =
      val lo = v.min
      val hi = v.max
      val pad = math.max((hi - lo) * padFrac, 1e-9)
      (lo - pad, hi + pad)
    val (xLo, xHi) = bounds(xs)
    val (yLo, yHi) = bounds(ys)

    val left = 52.0
    val right = 16.0
    val top = 14.0
    val bottom = 44.0
    val plotW = canvasWidth - left - right
    val plotH = canvasHeight - top - bottom
    def sx(v: Double): Double = left + (v - xLo) / (xHi - xLo) * plotW
    def sy(v: Double): Double = top + (yHi - v) / (yHi - yLo) * plotH

    // Recessive grid + tick labels.
    ctx.font = "11px system-ui, sans-serif"
    ctx.lineWidth = 1.0
    val xStep = niceStep((xHi - xLo) / 5.0)
    val yStep = niceStep((yHi - yLo) / 5.0)
    def ticks(lo: Double, hi: Double, step: Double): Seq[Double] =
      val first = math.ceil(lo / step) * step
      Iterator.iterate(first)(_ + step).takeWhile(_ <= hi + 1e-9).toSeq
    ctx.strokeStyle = th.grid
    ctx.fillStyle = th.textSecondary
    ctx.textAlign = "center"
    ctx.textBaseline = "top"
    for t <- ticks(xLo, xHi, xStep) do
      ctx.beginPath()
      ctx.moveTo(sx(t), top)
      ctx.lineTo(sx(t), top + plotH)
      ctx.stroke()
      ctx.fillText(formatTick(t, xStep), sx(t), top + plotH + 6)
    ctx.textAlign = "right"
    ctx.textBaseline = "middle"
    for t <- ticks(yLo, yHi, yStep) do
      ctx.beginPath()
      ctx.moveTo(left, sy(t))
      ctx.lineTo(left + plotW, sy(t))
      ctx.stroke()
      ctx.fillText(formatTick(t, yStep), left - 8, sy(t))

    // Plot frame (axes).
    ctx.strokeStyle = th.textSecondary
    ctx.strokeRect(left, top, plotW, plotH)

    // Axis titles.
    ctx.fillStyle = th.textPrimary
    ctx.font = "12px system-ui, sans-serif"
    ctx.textAlign = "center"
    ctx.textBaseline = "alphabetic"
    ctx.fillText("PC1", left + plotW / 2, canvasHeight - 10)
    ctx.save()
    ctx.translate(14, top + plotH / 2)
    ctx.rotate(-math.Pi / 2)
    ctx.fillText("PC2", 0, 0)
    ctx.restore()

    // Points: filled circles with a thin surface ring so overlaps stay legible.
    var i = 0
    while i < n do
      ctx.beginPath()
      ctx.arc(sx(xs(i)), sy(ys(i)), 4.0, 0, 2 * math.Pi)
      ctx.fillStyle = th.series(data.labels(i))
      ctx.fill()
      ctx.lineWidth = 1.0
      ctx.strokeStyle = th.surface
      ctx.stroke()
      i += 1

    // Direct cluster labels at the projected centroids (identity is never
    // carried by color alone).
    ctx.font = "600 12px system-ui, sans-serif"
    ctx.textAlign = "center"
    ctx.textBaseline = "middle"
    var c = 0
    while c < clusterNames.length do
      var cx = 0.0
      var cy = 0.0
      var count = 0
      var j = 0
      while j < n do
        if data.labels(j) == c then
          cx += xs(j)
          cy += ys(j)
          count += 1
        j += 1
      val label = clusterNames(c).takeRight(1) // "A" / "B" / "C"
      val px = sx(cx / count)
      val py = sy(cy / count)
      ctx.fillStyle = th.surface
      ctx.beginPath()
      ctx.arc(px, py, 9.0, 0, 2 * math.Pi)
      ctx.fill()
      ctx.strokeStyle = th.series(c)
      ctx.lineWidth = 1.5
      ctx.stroke()
      ctx.fillStyle = th.textPrimary
      ctx.fillText(label, px, py + 0.5)
      c += 1

  private def formatTick(v: Double, step: Double): String =
    if step >= 1.0 then math.round(v).toString
    else
      val decimals = math.max(0, -math.floor(math.log10(step)).toInt)
      val pow = math.pow(10.0, decimals.toDouble)
      val r = math.round(v * pow) / pow
      r.toString

  // ---------------------------------------------------------------------------
  // Page wiring.
  // ---------------------------------------------------------------------------

  private def pct(v: Double): String = f"${v * 100}%.1f%%"

  private var seed: Long = 42L
  private var lastData: Option[Dataset] = None
  private var lastResult: Option[PcaResult] = None

  private def setup(): Unit =
    val app = document.getElementById("app").asInstanceOf[html.Div]

    val controls = document.createElement("div").asInstanceOf[html.Div]
    controls.className = "controls"
    val button = document.createElement("button").asInstanceOf[html.Button]
    button.textContent = "Re-sample"
    button.title = "Regenerate the dataset and re-run the full PCA pipeline"
    val timing = document.createElement("span").asInstanceOf[html.Span]
    timing.className = "timing"
    controls.appendChild(button)
    controls.appendChild(timing)

    val canvas = document.createElement("canvas").asInstanceOf[html.Canvas]
    canvas.style.width = s"${canvasWidth}px"
    canvas.style.height = s"${canvasHeight}px"

    val legend = document.createElement("div").asInstanceOf[html.Div]
    legend.className = "legend"

    val readout = document.createElement("div").asInstanceOf[html.Div]
    readout.className = "readout"

    app.appendChild(controls)
    app.appendChild(canvas)
    app.appendChild(legend)
    app.appendChild(readout)

    def renderLegend(): Unit =
      legend.innerHTML = ""
      val th = theme
      clusterNames.zipWithIndex.foreach { case (name, idx) =>
        val item = document.createElement("span").asInstanceOf[html.Span]
        item.className = "legend-item"
        val dot = document.createElement("span").asInstanceOf[html.Span]
        dot.className = "legend-dot"
        dot.style.backgroundColor = th.series(idx)
        item.appendChild(dot)
        item.appendChild(document.createTextNode(name))
        legend.appendChild(item)
      }

    def recompute(): Unit =
      val data = generate(seed)
      runPca(data) match
        case Left(err) =>
          readout.textContent = s"PCA failed: $err"
          timing.textContent = ""
        case Right(result) =>
          lastData = Some(data)
          lastResult = Some(result)
          render(canvas, data, result)
          renderLegend()
          timing.textContent =
            f"center → scatter → eigSymmetric → project: ${result.elapsedMs}%.1f ms (n=${data.n}, d=$dims, seed=$seed)"
          readout.textContent =
            s"Variance explained — PC1: ${pct(result.varExplained1)}, PC2: ${pct(result.varExplained2)}, " +
              s"top-2 total: ${pct(result.varExplained1 + result.varExplained2)}"

    def redrawOnly(): Unit =
      for
        data <- lastData
        result <- lastResult
      do
        render(canvas, data, result)
        renderLegend()

    button.addEventListener("click", (_: dom.MouseEvent) => { seed += 1; recompute() })
    darkQuery.addEventListener("change", (_: dom.Event) => redrawOnly())

    recompute()

  def main(args: Array[String]): Unit =
    if document.readyState == dom.DocumentReadyState.loading then
      document.addEventListener("DOMContentLoaded", (_: dom.Event) => setup())
    else setup()
