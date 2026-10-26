import hashlib
import importlib.util
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
import zipfile

spec = importlib.util.spec_from_file_location('verify_bundle', Path(__file__).with_name('verify_bundle.py'))
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def jar(entry):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w') as z:
        z.writestr(entry, 'fixture')
        z.writestr('gale/Fixture.sjsir', 'fixture')
    return buf.getvalue()


def fixture():
    files = {}
    for a in v.ARTIFACTS:
        pom = f'''<project xmlns="http://maven.apache.org/POM/4.0.0">
        <modelVersion>4.0.0</modelVersion><groupId>{v.GROUP}</groupId><artifactId>{a}</artifactId>
        <version>0.1.0-M1</version><name>Gale</name><description>Fixture</description><url>https://example.com</url>
        <licenses><license><name>Apache-2.0</name><url>https://example.com/license</url></license></licenses>
        <developers><developer><id>dev</id><name>Developer</name></developer></developers>
        <scm><url>https://example.com</url><connection>scm:git:https://example.com</connection></scm>
        </project>'''.encode()
        base = v.PREFIX + a + '/0.1.0-M1/' + a + '-0.1.0-M1'
        for ext, data in [('.pom', pom), ('.jar', jar('gale/Fixture.class')),
                          ('-sources.jar', jar('Fixture.scala')), ('-javadoc.jar', jar('index.html'))]:
            files[base + ext] = data
    return files


def bundle(files, path):
    with zipfile.ZipFile(path, 'w') as z:
        for name, data in files.items():
            z.writestr(name, data)
            if not name.endswith('.asc'):
                for algo in ['md5', 'sha1']:
                    z.writestr(name + '.' + algo, hashlib.new(algo, data).hexdigest())


class BundleTests(unittest.TestCase):
    def check(self, mutate=None, valid=False):
        files = fixture()
        if mutate:
            mutate(files)
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / 'bundle.zip'
            bundle(files, path)
            if valid:
                self.assertEqual(v.verify(path), '0.1.0-M1')
            else:
                with self.assertRaises((ValueError, zipfile.BadZipFile, v.ET.ParseError)):
                    v.verify(path)

    def test_valid_unsigned(self):
        self.check(valid=True)

    def test_pom_metadata_and_coordinate_refusals(self):
        for old, new in [(b'<name>Gale</name>', b''), (b'<groupId>io.github.canardlapin</groupId>', b'<groupId>wrong</groupId>'),
                         (b'<version>0.1.0-M1</version>', b'<version>0.1.0-SNAPSHOT</version>')]:
            def mutate(files):
                name = next(n for n in files if n.endswith('.pom'))
                files[name] = files[name].replace(old, new)
            self.check(mutate)

    def test_snapshot_dependency(self):
        def mutate(files):
            name = next(n for n in files if n.endswith('.pom'))
            files[name] = files[name].replace(b'</project>', b'<dependencies><dependency><groupId>external</groupId>'
                                            b'<artifactId>bad</artifactId><version>1-SNAPSHOT</version></dependency></dependencies></project>')
        self.check(mutate)

    def test_invalid_archives_and_empty_payloads(self):
        for data in [b'not a jar', jar('irrelevant.txt')]:
            self.check(lambda files: files.update({next(n for n in files if n.endswith('-sources.jar')): data}))

    def test_missing_artifact_and_unexpected_file(self):
        self.check(lambda files: files.pop(next(iter(files))))
        self.check(lambda files: files.update({'evil/extra.jar': b'x'}))
        self.check(lambda files: files.update({'../outside': b'x'}))

    def test_checksums(self):
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / 'bundle.zip'
            bundle(fixture(), path)
            with zipfile.ZipFile(path) as z:
                files = {n: z.read(n) for n in z.namelist()}
            name = next(n for n in files if n.endswith('.sha1'))
            files[name] = b'0' * 40
            with zipfile.ZipFile(path, 'w') as z:
                for n, data in files.items():
                    z.writestr(n, data)
            with self.assertRaisesRegex(ValueError, 'incorrect sha1'):
                v.verify(path)

    def test_duplicate_paths(self):
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / 'bundle.zip'
            bundle(fixture(), path)
            with zipfile.ZipFile(path, 'a') as z:
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    z.writestr(z.namelist()[0], b'duplicate')
            with self.assertRaisesRegex(ValueError, 'duplicate archive paths'):
                v.verify(path)

    def test_signatures_and_expected_identity(self):
        with tempfile.TemporaryDirectory() as t:
            home = Path(t)
            def gpg(*args):
                return subprocess.run(['gpg', '--batch', '--homedir', str(home), *args], check=True, capture_output=True)
            gpg('--pinentry-mode', 'loopback', '--passphrase', '', '--quick-generate-key',
                'Gale disposable verifier test', 'ed25519', 'sign', '1d')
            fingerprint = next(line.split(':')[9] for line in gpg('--with-colons', '--list-keys').stdout.decode().splitlines()
                               if line.startswith('fpr:'))
            public = home / 'public.asc'
            public.write_bytes(gpg('--armor', '--export', fingerprint).stdout)
            files = fixture()
            for name, data in list(files.items()):
                payload = home / 'payload'
                payload.write_bytes(data)
                files[name + '.asc'] = gpg('--armor', '--detach-sign', '--output', '-', str(payload)).stdout
            path = home / 'bundle.zip'
            bundle(files, path)
            self.assertEqual(v.verify(path, True, public, fingerprint), '0.1.0-M1')
            with self.assertRaises(ValueError):
                v.verify(path, True, public, '0' * 40)
            files[next(n for n in files if n.endswith('.asc'))] = b'not a signature'
            bundle(files, path)
            with self.assertRaises(ValueError):
                v.verify(path, True, public, fingerprint)
            with self.assertRaises(ValueError):
                v.verify(path, True)


if __name__ == '__main__':
    unittest.main()
