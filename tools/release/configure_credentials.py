#!/usr/bin/env python3
"""Interactive owner setup; never prints secrets or writes exported keys to disk."""
import argparse
import base64
import getpass
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
REPO = 'canardlapin/gale'
UID = 'Gale releases <307091466+canardlapin@users.noreply.github.com>'
TOKEN_ENV = ('GH_TOKEN', 'GITHUB_TOKEN', 'GH_ENTERPRISE_TOKEN', 'GITHUB_ENTERPRISE_TOKEN')


def run(args, data=None, allow_failure=False):
    env = {k: v for k, v in os.environ.items() if k not in TOKEN_ENV}
    result = subprocess.run(args, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT, env=env)
    if result.returncode and not allow_failure:
        # In particular, never include private export output or secret-set input.
        raise RuntimeError(f'{args[0]} {args[1]} failed (exit {result.returncode}); no secret output displayed')
    return result.stdout if result.returncode == 0 else b''


def account_check():
    if run(['gh', 'api', 'user', '--jq', '.login']).decode().strip() != 'canardlapin':
        raise RuntimeError('Use the repo-local canardlapin GitHub CLI account; no credentials were changed')


def fingerprints(selector):
    listing = run(['gpg', '--batch', '--with-colons', '--list-secret-keys', selector], allow_failure=True).decode()
    found, primary = [], False
    for line in listing.splitlines():
        fields = line.split(':')
        if fields[0] == 'sec':
            primary = True
        elif fields[0] == 'ssb':
            primary = False
        elif fields[0] == 'fpr' and primary:
            found.append(fields[9])
            primary = False
    return found


def set_secret(name, value):
    run(['gh', 'secret', 'set', name, '--repo', REPO], value)
    print(f'Configured {name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Check public configuration names only; no changes')
    parser.add_argument('--fingerprint', help='Reuse this existing signing key instead of the dedicated Gale key')
    args = parser.parse_args()
    account_check()
    if args.check:
        for kind in ('secret', 'variable'):
            names = json.loads(run(['gh', kind, 'list', '--repo', REPO, '--json', 'name']))
            print(kind + ' names: ' + ', '.join(item['name'] for item in names))
        print('Dedicated local key fingerprints: ' + ', '.join(fingerprints(UID)))
        return
    if not sys.stdin.isatty():
        raise RuntimeError('Run this command directly in Terminal; credentials require hidden interactive input')
    if args.fingerprint and not re.fullmatch(r'[0-9A-Fa-f]{40}|[0-9A-Fa-f]{64}', args.fingerprint):
        raise RuntimeError('Supply a full primary-key fingerprint')
    print('This configures GitHub release secrets, creates/reuses an encrypted local signing key,')
    print('and sends its PUBLIC key to keyserver.ubuntu.com. It does not publish a Gale release.')
    print('First sign into https://central.sonatype.com as GitHub user canardlapin, verify')
    print('io.github.canardlapin, and create a token at https://central.sonatype.com/usertoken.')
    keys = fingerprints(args.fingerprint or UID)
    if len(keys) > 1:
        raise RuntimeError('Multiple matching keys; rerun with --fingerprint to select one')
    if args.fingerprint and not keys:
        raise RuntimeError('The selected secret key is not available locally')
    passphrase = getpass.getpass('Signing-key passphrase (save it in your password manager): ')
    if not passphrase or '\n' in passphrase or '\r' in passphrase:
        raise RuntimeError('A nonempty single-line signing passphrase is required')
    if not keys and passphrase != getpass.getpass('Confirm new signing passphrase: '):
        raise RuntimeError('Passphrases do not match')
    username = getpass.getpass('Central token username (hidden): ').strip()
    password = getpass.getpass('Central token password (hidden): ').strip()
    if not username or not password:
        raise RuntimeError('Both Central token fields are required')
    phrase = (passphrase + '\n').encode()
    if not keys:
        print('Generating the dedicated RSA-4096 signing key (expires in two years)...')
        run(['gpg', '--batch', '--pinentry-mode', 'loopback', '--passphrase-fd', '0',
             '--quick-generate-key', UID, 'rsa4096', 'sign', '2y'], phrase)
        keys = fingerprints(UID)
    if len(keys) != 1:
        raise RuntimeError('Cannot identify the dedicated signing key')
    fingerprint = keys[0]
    private = run(['gpg', '--batch', '--pinentry-mode', 'loopback', '--passphrase-fd', '0',
                   '--armor', '--export-secret-keys', fingerprint], phrase)
    if not private.startswith(b'-----BEGIN PGP PRIVATE KEY BLOCK-----'):
        raise RuntimeError('Private-key export failed')
    run(['gpg', '--batch', '--keyserver', 'hkps://keyserver.ubuntu.com', '--send-keys', fingerprint])
    set_secret('PGP_SECRET', base64.b64encode(private))
    set_secret('PGP_PASSPHRASE', passphrase.encode())
    set_secret('SONATYPE_USERNAME', username.encode())
    set_secret('SONATYPE_PASSWORD', password.encode())
    run(['gh', 'variable', 'set', 'RELEASE_PGP_FINGERPRINT', '--repo', REPO, '--body', fingerprint])
    print('Configured public signing fingerprint: ' + fingerprint)
    print('The encrypted private key stays in your GPG keyring; exported key material was piped directly to GitHub.')
    print('Keep the passphrase and Portal token in your password manager. Back up your encrypted GPG keyring and revocation certificate.')
    print('Configuration is complete; real signed Portal validation is still required before release.')


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, KeyboardInterrupt, EOFError) as error:
        sys.exit(str(error) or 'Setup cancelled; any previously completed steps remain in place.')
