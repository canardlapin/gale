// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import { repo, sidebar } from './docs-manifest.mjs';

// GitHub Pages serves the project site under /gale. SITE_BASE overrides it for
// another host; scripts/sync-docs.mjs reads the same variable when it rewrites links.
const base = process.env.SITE_BASE ?? '/gale';

export default defineConfig({
	site: 'https://canardlapin.github.io',
	base,
	trailingSlash: 'always',
	integrations: [
		starlight({
			title: 'Gale',
			description:
				'Linear algebra for Scala 3 on the JVM and Scala.js: dense and sparse matrices, solvers, spectral methods, and optimization, with typed failures and convergence diagnostics.',
			logo: { src: './src/assets/mark.svg', alt: '' },
			favicon: '/favicon.svg',
			social: [{ icon: 'github', label: 'GitHub', href: repo }],
			customCss: [
				'@fontsource-variable/newsreader/opsz.css',
				'@fontsource-variable/newsreader/opsz-italic.css',
				'@fontsource-variable/inter',
				'@fontsource-variable/jetbrains-mono',
				'./src/styles/theme.css',
			],
			expressiveCode: {
				themes: ['github-light', 'github-dark-dimmed'],
				styleOverrides: {
					borderRadius: '0.4rem',
					codeFontFamily: "'JetBrains Mono Variable', ui-monospace, monospace",
					codeFontSize: '0.84rem',
				},
			},
			lastUpdated: true,
			sidebar: sidebar(),
		}),
	],
});
