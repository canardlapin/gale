import importlib.util
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('configure', Path(__file__).with_name('configure_credentials.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CredentialTests(unittest.TestCase):
    def test_secret_uses_stdin_and_token_overrides_are_removed(self):
        with patch.dict(c.os.environ, {'GH_TOKEN': 'synthetic-override'}), patch.object(c.subprocess, 'run') as proc:
            proc.return_value = subprocess.CompletedProcess([], 0, b'', b'')
            c.set_secret('PGP_SECRET', b'synthetic-secret')
            args, kwargs = proc.call_args
            self.assertNotIn('synthetic-secret', str(args))
            self.assertEqual(kwargs['input'], b'synthetic-secret')
            self.assertNotIn('GH_TOKEN', kwargs['env'])

    def test_failed_export_does_not_include_secret_in_error(self):
        with patch.object(c.subprocess, 'run') as proc:
            proc.return_value = subprocess.CompletedProcess([], 1, b'private-key', b'sensitive-error')
            with self.assertRaises(RuntimeError) as error:
                c.run(['gpg', '--export-secret-keys'], b'passphrase')
            self.assertNotIn('private-key', str(error.exception))
            self.assertNotIn('sensitive-error', str(error.exception))
            self.assertNotIn('passphrase', str(error.exception))

    def test_wrong_github_account_refuses(self):
        with patch.object(c, 'run', return_value=b'wrong-account\n'):
            with self.assertRaisesRegex(RuntimeError, 'no credentials were changed'):
                c.account_check()


if __name__ == '__main__':
    unittest.main()
