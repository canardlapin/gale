"""Exercise the actual publication workflow without network or credentials."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import textwrap
import unittest
from unittest.mock import patch
import urllib.request

WORKFLOW = Path(__file__).resolve().parents[2] / '.github/workflows/publish-central.yml'
CODE = textwrap.dedent(WORKFLOW.read_text().split('        run: |\n', 1)[1].split('      - uses:', 1)[0])
DEPLOYMENT = '00000000-0000-0000-0000-000000000001'
PURLS = ['pkg:maven/io.github.canardlapin/' + artifact + '@0.1.0-M1'
         for artifact in ('gale-core_3', 'gale-laws_3', 'gale-core_sjs1_3', 'gale-laws_sjs1_3')]


class Response:
    def __init__(self, data):
        self.data = data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def read(self):
        return self.data


class PublicationWorkflowTest(unittest.TestCase):
    def run_states(self, states, publishes, error=None):
        requests = []
        sequence = iter(states)

        def request(req, timeout):
            requests.append(req.full_url)
            if '/deployment/' in req.full_url:
                return Response(b'')
            state, purls = next(sequence)
            return Response(json.dumps(dict(deploymentId=DEPLOYMENT, deploymentState=state, purls=purls)).encode())

        with tempfile.TemporaryDirectory() as temporary:
            old = os.getcwd()
            os.chdir(temporary)
            try:
                with patch.dict(os.environ, dict(DEPLOYMENT_ID=DEPLOYMENT, SONATYPE_USERNAME='synthetic',
                                                SONATYPE_PASSWORD='synthetic', GITHUB_STEP_SUMMARY=temporary + '/summary')), \
                     patch.object(urllib.request, 'urlopen', request), patch('time.sleep'), \
                     contextlib.redirect_stdout(io.StringIO()):
                    if error:
                        with self.assertRaises(error):
                            exec(compile(CODE, str(WORKFLOW), 'exec'), {})
                    else:
                        exec(compile(CODE, str(WORKFLOW), 'exec'), {})
            finally:
                os.chdir(old)
        self.assertEqual(sum('/deployment/' in url for url in requests), publishes)

    def test_central_can_omit_purls_during_publication(self):
        self.run_states([('VALIDATED', PURLS), ('PUBLISHING', []), ('PUBLISHED', PURLS)], 1)

    def test_resume_does_not_publish_again(self):
        self.run_states([('PUBLISHING', []), ('PUBLISHED', PURLS)], 0)

    def test_missing_prepublication_coordinates_are_rejected(self):
        self.run_states([('VALIDATED', [])], 0, AssertionError)

    def test_wrong_coordinates_are_rejected_even_during_publication(self):
        self.run_states([('PUBLISHING', PURLS[:1])], 0, AssertionError)

    def test_completion_requires_exact_coordinates(self):
        self.run_states([('PUBLISHED', [])], 0, AssertionError)

    def test_failed_validation_cannot_publish(self):
        self.run_states([('FAILED', PURLS)], 0, RuntimeError)


if __name__ == '__main__':
    unittest.main()
