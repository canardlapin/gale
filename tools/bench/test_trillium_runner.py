#!/usr/bin/env python3
"""Exercise runner argument/lifecycle contracts without Slurm, Java or Linux."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


RUNNER = Path(__file__).resolve().parents[2] / "docs/verification/w21-simd-spike/x86-avx512/run.sh"
SHA = "41cbf87505af14c3af36eab1563b2ba42cea7a61"


class TrilliumRunnerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bundle = self.root / "bundle"
        self.bundle.mkdir()
        (self.bundle / "SHA").write_text(SHA + "\n")
        (self.bundle / "jmh.cp").write_text("fake-classpath\n")
        self.out = self.root / "receipt"
        self.out.mkdir()
        self.calls = self.root / "calls.jsonl"
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.script("module", '#!/bin/bash\nexit "${MODULE_EXIT:-0}"\n')
        self.script("cat", f'''#!{sys.executable}
import pathlib, sys
for path in sys.argv[1:]:
    if path.startswith('/sys/devices/system/cpu/'):
        print('0-7')
    else:
        print(pathlib.Path(path).read_text(), end='')
''')
        self.script("ls", '#!/bin/bash\necho /sys/devices/system/cpu/cpu0/node0\n')
        for command in ("lscpu", "ldconfig", "uptime"):
            self.script(command, '#!/bin/bash\necho synthetic-host\n')
        self.script("numactl", '''#!/bin/bash
if [[ "$1" == "-H" ]]; then exit 0; fi
shift 2
exec "$@"
''')
        self.script("java", f'''#!{sys.executable}
import json, os, sys
if '-version' in sys.argv:
    print('synthetic java 25', file=sys.stderr)
    sys.exit(0)
with open(os.environ['CALLS'], 'a') as f:
    f.write(json.dumps({{'args': sys.argv[1:], 'inherited': [os.environ.get(k) for k in
        ('JAVA_TOOL_OPTIONS', 'JDK_JAVA_OPTIONS', '_JAVA_OPTIONS')]}}) + '\\n')
sys.exit(int(os.environ.get('JAVA_EXIT', '0')))
''')
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ["PATH"],
                        GALE_BENCH_ROOT=str(self.root), CALLS=str(self.calls),
                        JAVA_TOOL_OPTIONS="-Xmx2g", JDK_JAVA_OPTIONS="-Xmx2g",
                        _JAVA_OPTIONS="-Xmx2g", LD_LIBRARY_PATH="/fake/library/path")

    def script(self, name, content):
        path = self.bin / name
        path.write_text(content)
        path.chmod(0o755)

    def run_mode(self, mode):
        return subprocess.run(["bash", str(RUNNER), mode, str(self.out)],
                              env=self.env, text=True, capture_output=True)

    def recorded(self):
        return [json.loads(line) for line in self.calls.read_text().splitlines()] if self.calls.exists() else []

    def test_both_lanes_bind_source_and_reset_both_sidecars(self):
        for mode, option in (("laneA", "-jvmArgs"), ("laneB", "-jvmArgsAppend")):
            with self.subTest(mode=mode):
                for name in ("netlib.jsonl", "netlib.jsonl.validation.jsonl", "result.json", "keep.txt"):
                    (self.out / name).write_text("stale")
                result = self.run_mode(mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                call = self.recorded()[-1]
                args = call["args"]
                self.assertEqual(args[args.index("-foe") + 1], "true")
                fork = args[args.index(option) + 1]
                self.assertIn("-Dgale.bench.sourceSha=" + SHA, fork)
                self.assertIn("-Dgale.bench.netlibSidecar=" + str(self.out / "netlib.jsonl"), fork)
                self.assertIn("-Dgale.bench.lane=" + mode[-1], fork)
                self.assertEqual(call["inherited"], [None, None, None])
                for name in ("netlib.jsonl", "netlib.jsonl.validation.jsonl", "result.json"):
                    self.assertFalse((self.out / name).exists())
                self.assertEqual((self.out / "keep.txt").read_text(), "stale")
                self.assertEqual((self.out / "rc.txt").read_text(), "rc=0\n")

    def test_env_probe_preserves_receipts(self):
        for name in ("netlib.jsonl", "netlib.jsonl.validation.jsonl", "result.json"):
            (self.out / name).write_text("preserve")
        result = self.run_mode("env")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.recorded(), [])
        self.assertEqual((self.out / "netlib.jsonl.validation.jsonl").read_text(), "preserve")

    def test_bad_sha_fails_before_jmh(self):
        (self.bundle / "SHA").write_text("short-sha\n")
        result = self.run_mode("laneA")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(self.recorded(), [])
        self.assertEqual((self.out / "rc.txt").read_text(), "rc=2\n")

    def test_missing_module_fails_before_jmh(self):
        self.env["MODULE_EXIT"] = "8"
        result = self.run_mode("laneB")
        self.assertEqual(result.returncode, 8)
        self.assertEqual(self.recorded(), [])
        self.assertEqual((self.out / "rc.txt").read_text(), "rc=8\n")

    def test_gate_failure_is_not_hidden_by_second_command(self):
        self.env["JAVA_EXIT"] = "7"
        result = self.run_mode("gate")
        self.assertEqual(result.returncode, 7, result.stderr)
        self.assertEqual(len(self.recorded()), 1)
        self.assertIn("-Dgale.bench.sourceSha=" + SHA, self.recorded()[0]["args"])
        self.assertEqual((self.out / "rc.txt").read_text(), "rc=7\n")


if __name__ == "__main__":
    unittest.main()
