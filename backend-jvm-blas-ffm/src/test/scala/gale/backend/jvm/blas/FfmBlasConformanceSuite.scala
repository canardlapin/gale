package gale.backend.jvm.blas

import gale.backend.Backend
import gale.laws.BackendConformanceSuite

class FfmBlasConformanceSuite extends BackendConformanceSuite:
  private lazy val loaded = FfmBlasBackend.load().fold(throw _, identity)
  def backend: Backend = loaded

  override def afterAll(): Unit =
    if loaded != null then loaded.close()
