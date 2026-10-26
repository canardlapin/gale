#!/usr/bin/env python3
"""Validate the four-artifact M1 bundle before Central Portal upload."""
import argparse
import hashlib
import io
from pathlib import Path, PurePosixPath
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zipfile

ARTIFACTS = {'gale-core_3', 'gale-core_sjs1_3', 'gale-laws_3', 'gale-laws_sjs1_3'}
GROUP = 'io.github.canardlapin'
PREFIX = 'io/github/canardlapin/'
VERSION = r'0\.1\.\d+(?:-(?:M|RC)[1-9]\d*)?'
NS = {'m': 'http://maven.apache.org/POM/4.0.0'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_pom(data, artifact, version):
    root = ET.fromstring(data)
    def value(path):
        return (root.findtext(path, namespaces=NS) or '').strip()
    for field, expected in [('groupId', GROUP), ('artifactId', artifact), ('version', version)]:
        require(value('m:' + field) == expected, f'{artifact}: POM {field} does not match its path')
    for path in ['m:modelVersion', 'm:name', 'm:description', 'm:url',
                 'm:licenses/m:license/m:name', 'm:licenses/m:license/m:url',
                 'm:developers/m:developer/m:id', 'm:developers/m:developer/m:name',
                 'm:scm/m:url', 'm:scm/m:connection']:
        require(value(path), f'{artifact}: missing POM metadata {path}')
    for dep in root.findall('.//m:dependency', NS):
        group, name, revision = [(dep.findtext('m:' + k, namespaces=NS) or '').strip()
                                 for k in ['groupId', 'artifactId', 'version']]
        require(group and name and revision, f'{artifact}: incomplete dependency')
        require('SNAPSHOT' not in revision.upper() and '${' not in revision and
                not any(c in revision for c in '[](),+') and revision not in {'LATEST', 'RELEASE'},
                f'{artifact}: non-immutable dependency {name}:{revision}')
        if group == GROUP:
            require(name in ARTIFACTS and revision == version,
                    f'{artifact}: dependency outside admitted candidate: {name}:{revision}')
    require(not root.findall('m:repositories', NS), f'{artifact}: unexpected repositories')


def validate_jar(data, suffix, artifact):
    with zipfile.ZipFile(io.BytesIO(data)) as jar:
        require(jar.testzip() is None, f'{artifact}: corrupt {suffix}')
        names = jar.namelist()
        ending = '.scala' if suffix == '-sources.jar' else '.html' if suffix == '-javadoc.jar' else '.sjsir' if artifact.endswith('_sjs1_3') else '.class'
        require(any(n.endswith(ending) for n in names), f'{artifact}: empty or wrong {suffix} contents')
        if artifact.endswith('_sjs1_3') and suffix == '.jar':
            require(any(n.endswith('.sjsir') for n in names), f'{artifact}: missing Scala.js IR')


def verify(bundle, require_signatures=False, public_key=None, expected_fingerprint=None):
    if require_signatures:
        require(public_key and expected_fingerprint, 'signed verification requires a public key and expected fingerprint')
        require(re.fullmatch(r'[0-9A-Fa-f]{40}|[0-9A-Fa-f]{64}', expected_fingerprint), 'invalid signing fingerprint')
    with zipfile.ZipFile(bundle) as archive, tempfile.TemporaryDirectory() as scratch:
        paths = archive.namelist()
        require(len(paths) == len(set(paths)), 'duplicate archive paths')
        files = {n: archive.read(n) for n in paths if not n.endswith('/')}
        for n in paths:
            path = PurePosixPath(n)
            require(not path.is_absolute() and '..' not in path.parts and '\\' not in n, 'unsafe archive path')
        poms = [n for n in files if n.endswith('.pom')]
        require(len(poms) == 4, 'expected exactly four POMs')
        version = None
        seen = set()
        allowed = set()
        home = Path(scratch)
        if require_signatures:
            imported = subprocess.run(['gpg', '--batch', '--homedir', str(home), '--import', str(public_key)],
                                      capture_output=True)
            require(imported.returncode == 0, 'cannot import release public key')
        for pom in poms:
            require(pom.startswith(PREFIX), 'coordinate outside admitted group')
            parts = pom[len(PREFIX):].split('/')
            require(len(parts) == 3, 'invalid Maven layout')
            artifact, revision, filename = parts
            require(artifact in ARTIFACTS and artifact not in seen, 'unexpected or duplicate coordinate')
            require(re.fullmatch(VERSION, revision), 'invalid release version')
            require(version is None or revision == version, 'mixed versions')
            version = revision
            seen.add(artifact)
            base = PREFIX + artifact + '/' + revision + '/' + artifact + '-' + revision
            require(pom == base + '.pom', 'POM filename mismatch')
            validate_pom(files[pom], artifact, version)
            for suffix in ['.pom', '.jar', '-sources.jar', '-javadoc.jar']:
                name = base + suffix
                require(name in files, f'missing {name}')
                allowed.add(name)
                data = files[name]
                if suffix != '.pom':
                    validate_jar(data, suffix, artifact)
                for algorithm in ['md5', 'sha1']:
                    checksum = name + '.' + algorithm
                    require(checksum in files and files[checksum].decode().strip() == hashlib.new(algorithm, data).hexdigest(),
                            f'missing or incorrect {algorithm}: {name}')
                    allowed.add(checksum)
                signature = name + '.asc'
                if signature in files:
                    allowed.add(signature)
                    for algorithm in ['md5', 'sha1']:
                        checksum = signature + '.' + algorithm
                        if checksum in files:
                            require(files[checksum].decode().strip() == hashlib.new(algorithm, files[signature]).hexdigest(),
                                    f'incorrect signature checksum: {name}')
                            allowed.add(checksum)
                if require_signatures:
                    require(signature in files, f'missing signature: {name}')
                    payload, detached = home / 'payload', home / 'payload.asc'
                    payload.write_bytes(data)
                    detached.write_bytes(files[signature])
                    verified = subprocess.run(['gpg', '--batch', '--homedir', str(home), '--status-fd', '1',
                                               '--verify', str(detached), str(payload)], capture_output=True, text=True)
                    valid = [line.split()[2:] for line in verified.stdout.splitlines()
                             if line.startswith('[GNUPG:] VALIDSIG ')]
                    fingerprint = expected_fingerprint.upper()
                    require(verified.returncode == 0 and any(v[0] == fingerprint or v[-1] == fingerprint for v in valid),
                            f'invalid signature or unexpected signer: {name}')
        require(seen == ARTIFACTS, 'artifact set mismatch')
        require(set(files) == allowed, f'unexpected bundle files: {sorted(set(files) - allowed)[:3]}')
        return version


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-signatures', action='store_true')
    parser.add_argument('--public-key', type=Path)
    parser.add_argument('--expected-fingerprint')
    parser.add_argument('bundle', type=Path)
    args = parser.parse_args()
    try:
        version = verify(args.bundle, args.require_signatures, args.public_key, args.expected_fingerprint)
    except (ValueError, OSError, ET.ParseError, zipfile.BadZipFile) as error:
        parser.exit(1, f'Central bundle rejected: {error}\n')
    print(f'Central bundle verified: 4 admitted artifacts, version {version}; signatures ' +
          ('authenticated' if args.require_signatures else 'not authenticated (unsigned rehearsal)'))


if __name__ == '__main__':
    main()
