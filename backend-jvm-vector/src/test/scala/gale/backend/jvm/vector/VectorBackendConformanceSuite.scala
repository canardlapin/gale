package gale.backend.jvm.vector

import gale.backend.Backend
import gale.laws.BackendConformanceSuite

/** Runs the same reusable contract that every optional accelerator must pass. */
class VectorBackendConformanceSuite extends BackendConformanceSuite:
  val backend: Backend = VectorBackend
