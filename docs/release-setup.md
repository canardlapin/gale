# Set up Gale release credentials

Use the `canardlapin` GitHub identity and namespace `io.github.canardlapin`.
This setup stores credentials; it does not create a release tag, upload a bundle,
or publish Gale.

1. Sign into [Central Portal](https://central.sonatype.com) using GitHub as
   **canardlapin**. Check **View Namespaces** for `io.github.canardlapin`.
   GitHub sign-up normally provisions the matching namespace automatically;
   otherwise follow [Sonatype's namespace instructions](https://central.sonatype.org/register/namespace/).
2. Open [User Tokens](https://central.sonatype.com/usertoken), generate a token
   named `gale-github-actions`, and choose an expiration. Save both generated
   fields in your password manager before closing the modal. These are the
   **token username and token password**, not your GitHub login/password.
   See [Sonatype's token guide](https://central.sonatype.org/publish/generate-portal-token/).
3. In your own Terminal, run:

   ```sh
   cd /Users/bbuchsbaum/code/scala/gale
   python3 tools/release/configure_credentials.py
   ```

   Choose and save a signing passphrase, then enter the two Portal token fields
   at the hidden prompts. Do not paste credentials into chat or shell commands.

The helper checks that the repository's account-aware `gh` authenticates as
`canardlapin`. It creates a dedicated RSA-4096 signing key expiring in two years
(or reuses the matching dedicated key), publishes its public key to the
Sonatype-supported Ubuntu keyserver, and sets:

| GitHub setting | Kind |
|---|---|
| `PGP_SECRET` | Secret: base64-encoded, passphrase-protected private key |
| `PGP_PASSPHRASE` | Secret |
| `SONATYPE_USERNAME` | Secret: Portal token username |
| `SONATYPE_PASSWORD` | Secret: Portal token password |
| `RELEASE_PGP_FINGERPRINT` | Variable: public primary-key fingerprint |

Secret values travel through hidden input and subprocess pipes, not command-line
arguments, repository files, or printed logs. The encrypted private key remains
in your local GnuPG keyring. Back up that keyring and its revocation certificate
securely, and retain the passphrase and Portal token in your password manager.
If setup is interrupted, rerun it: the dedicated key is reused. Use
`--fingerprint FULL_PRIMARY_FINGERPRINT` only to deliberately select another
existing signing key. Publishing the public key follows
[Sonatype's signing requirements](https://central.sonatype.org/publish/requirements/gpg/).

Check names and public fingerprints without revealing values:

```sh
python3 tools/release/configure_credentials.py --check
```

Then qualify the exact candidate with CI and the unsigned Release workflow dry
run. A real signed deployment must reach **VALIDATED** in USER_MANAGED mode.
The workflow records the validation response and does not request publication.
Final Central publication remains a separate owner-authorized step.
