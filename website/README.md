# Gale documentation site

The [Starlight](https://starlight.astro.build) site for
<https://canardlapin.github.io/gale/>.

The guides are not written here. They are the mdoc output of
[`docs/user`](../docs/user), so every example on the site was compiled and run
against the current `gale-core`. `scripts/sync-docs.mjs` copies that output into
`src/content/docs/`, moves each page's H1 into frontmatter, rewrites links to site
routes (and links outside `docs/user` to GitHub), and drops the object-hash lines
mdoc prints for values without a useful `toString`. Only the landing page,
`src/content/docs/index.mdx`, is written in this directory.

| To change | Edit |
| --- | --- |
| A guide | its source in `docs/user/` |
| The sidebar, labels, or URLs | `docs-manifest.mjs` |
| The landing page | `src/content/docs/index.mdx` |
| The code shown on the landing page | the guide section named in `landingSnippets` in `scripts/sync-docs.mjs` |
| The landing figure | `scripts/draw-hero.mjs`, then run it |
| Colours and type | `src/styles/theme.css` |

The sync fails if a guide in `docs/user` is missing from `docs-manifest.mjs`, if a
relative link escapes the rewriter, or if a landing-page section has been renamed.

This directory is `website/`, not `site/`, because `site/` is the sbt project
that runs mdoc and the Laika guide.

## Build

Requires Node 22.23.4 or later.

```sh
sbt docs/mdoc            # in the repository root: compile and run every guide example
cd website
npm ci
npm run build            # refuses mdoc output older than its docs/user source
npm run preview
```

`npm run dev` syncs with `--allow-stale`, so you can work on the landing page or
the theme without recompiling. It warns about stale pages but does not stop.

`.github/workflows/docs-site.yml` runs the same steps on every pull request that
touches the guides, and deploys `dist/` to GitHub Pages from `main`.
