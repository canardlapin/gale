#!/usr/bin/env node
// Copy the compiled guides into the site.
//
// The source is mdoc's output (site/target/mdoc), not docs/user itself: mdoc has
// already compiled every example and appended the values it printed. This script only
// re-homes those pages. It moves each H1 into frontmatter, rewrites links to the
// site's routes, and drops the identity-hash lines that mdoc prints for values with no
// useful toString.
//
//   node scripts/sync-docs.mjs                # requires a current `sbt docs/mdoc`
//   node scripts/sync-docs.mjs --allow-stale  # draft against older mdoc output

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { allPages, branch, repo } from '../docs-manifest.mjs';

const siteRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const repoRoot = path.resolve(siteRoot, '..');
const docsSource = path.join(repoRoot, 'docs/user');
const compiled = path.join(repoRoot, 'site/target/mdoc');
const contentOut = path.join(siteRoot, 'src/content/docs');
const generatedOut = path.join(siteRoot, 'src/generated');
const base = (process.env.SITE_BASE ?? '/gale').replace(/\/$/, '');
const allowStale = process.argv.includes('--allow-stale');

const fail = (message) => {
	console.error(`sync-docs: ${message}`);
	process.exit(1);
};

if (!fs.existsSync(compiled)) {
	fail(`${path.relative(repoRoot, compiled)} does not exist; run \`sbt docs/mdoc\` in the repository root first`);
}

const pages = allPages();
const bySource = new Map(pages.map((page) => [page.source, page]));

// Every guide must have a place on the site, or it would silently go unpublished.
const listed = new Set(pages.map((page) => page.source));
const walk = (dir) =>
	fs.readdirSync(dir, { withFileTypes: true }).flatMap((entry) => {
		const full = path.join(dir, entry.name);
		return entry.isDirectory() ? walk(full) : entry.name.endsWith('.md') ? [path.relative(docsSource, full)] : [];
	});
const unlisted = walk(docsSource).filter((source) => !listed.has(source.split(path.sep).join('/')));
if (unlisted.length > 0) fail(`these guides are missing from docs-manifest.mjs:\n${unlisted.map((s) => `  docs/user/${s}`).join('\n')}`);

// A page whose source is newer than its mdoc output was edited after the last compile,
// so the site would publish examples that were never run.
const stale = pages.filter((page) => {
	const src = path.join(docsSource, page.source);
	const out = path.join(compiled, page.source);
	if (!fs.existsSync(src)) fail(`docs/user/${page.source} is in the page map but does not exist`);
	if (!fs.existsSync(out)) return true;
	return fs.statSync(src).mtimeMs > fs.statSync(out).mtimeMs;
});
if (stale.length > 0) {
	const list = stale.map((page) => `  docs/user/${page.source}`).join('\n');
	const message = `mdoc output is older than these sources; run \`sbt docs/mdoc\`:\n${list}`;
	if (allowStale) console.warn(`sync-docs: warning: ${message}`);
	else fail(message);
}

/** GitHub URL for a repository path outside the published page map. */
function githubUrl(repoPath) {
	const full = path.join(repoRoot, repoPath);
	if (!fs.existsSync(full)) fail(`a guide links to ${repoPath}, which does not exist`);
	const kind = fs.statSync(full).isDirectory() ? 'tree' : 'blob';
	return `${repo}/${kind}/${branch}/${repoPath}`;
}

/** Rewrite one Markdown link target found in docs/user/<source>. */
function rewriteTarget(target, source) {
	if (/^[a-z][a-z0-9+.-]*:/i.test(target) || target.startsWith('#')) return target;
	const [rawPath, anchor] = target.split('#', 2);
	const hash = anchor === undefined ? '' : `#${anchor}`;
	const fromDocs = path.posix.normalize(path.posix.join(path.posix.dirname(source), rawPath));

	const page = bySource.get(fromDocs);
	if (page) return `${base}/${page.slug}/${hash}`;

	const repoPath = path.posix.normalize(path.posix.join('docs/user', fromDocs));
	if (repoPath.startsWith('..')) fail(`docs/user/${source} links outside the repository: ${target}`);
	return githubUrl(repoPath) + hash;
}

/** Split Markdown into prose and fenced-code segments. */
function segments(body) {
	const out = [];
	let lines = [];
	let fence = null;
	const flush = (kind) => {
		if (lines.length > 0) out.push({ kind, text: lines.join('\n') });
		lines = [];
	};
	for (const line of body.split('\n')) {
		const open = line.match(/^ {0,3}(`{3,}|~{3,})/);
		if (fence === null && open) {
			flush('prose');
			fence = open[1];
			lines.push(line);
		} else if (fence !== null) {
			lines.push(line);
			// A fence closes on a bare run of the same character, at least as long as the opener.
			const close = line.match(/^ {0,3}(`{3,}|~{3,})\s*$/);
			if (close && close[1][0] === fence[0] && close[1].length >= fence.length) {
				fence = null;
				flush('code');
			}
		} else {
			lines.push(line);
		}
	}
	flush(fence === null ? 'prose' : 'code');
	return out;
}

/** Rewrite inline links outside code. Link text may span lines but not paragraphs. */
function rewriteLinks(prose, source) {
	// An inline code span is matched first and returned unchanged, so a `[a](b)` written
	// as code stays literal while [`Code`](b.md) is still a link.
	const token = /(?<!`)(`+)(?!`)[\s\S]*?(?<!`)\1(?!`)|(!?\[(?:[^\]`]|`[^`]*`)*\])\(([^)\s]+)\)/g;
	return prose
		.split(/(\n[ \t]*\n)/)
		.map((part, i) =>
			i % 2 === 1
				? part
				: part.replace(token, (match, _code, text, target) =>
						text === undefined ? match : `${text}(${rewriteTarget(target, source)})`,
					),
		)
		.join('');
}

// mdoc prints `// name: DMat = gale.linalg.DMat@2c5c2503` for a value whose type
// does not override toString. The hash changes on every run and tells the reader nothing.
const identityHash = /^\/\/ [\w`$]+: [\w.\[\], ]+ = [\w.$]+@[0-9a-f]+$/;
const tidyCode = (code) =>
	code
		.split('\n')
		.filter((line) => !identityHash.test(line))
		.join('\n');

/** First prose paragraph, flattened to plain text, for the page description. */
function describe(body) {
	const paragraph = body
		.split(/\n\s*\n/)
		.map((block) => block.trim())
		.find((block) => block && !/^(#|```|~~~|\||<|>|-\s|\*\s|\d+\.\s|!\[)/.test(block));
	if (!paragraph) return undefined;
	const text = paragraph
		.replace(/!?\[([^\]]*)\]\([^)]*\)/g, '$1')
		.replace(/`|\*\*?(?=\S)|(?<=\S)\*\*?/g, '')
		.replace(/\s+/g, ' ')
		.trim();
	if (text.length <= 200) return text;
	const cut = text.slice(0, 200);
	return `${cut.slice(0, cut.lastIndexOf(' '))}…`;
}

const yaml = (value) => JSON.stringify(value);

// Remove the previous sync so a page dropped from the map does not linger.
for (const entry of fs.readdirSync(contentOut)) {
	if (entry === 'index.mdx') continue;
	fs.rmSync(path.join(contentOut, entry), { recursive: true, force: true });
}

/** Compiled code shown on the landing page, keyed by name. */
const snippets = {};

/** The sections whose code the landing page shows, exactly as mdoc compiled it. */
const landingSnippets = [
	{ name: 'first-program', source: 'getting-started.md', headings: ['Make a matrix and use it'] },
	{ name: 'solve', source: 'getting-started.md', headings: ['Solve a system', 'Keep numerical failure explicit'] },
]

/**
 * A landing-page excerpt keeps the code and the final printed result. It drops the
 * guide's assertions and any printed value that still shows an object hash, so the
 * reader sees the program and the answer rather than the test scaffolding.
 */
function forLanding(code) {
	const out = [];
	let block = [];
	const flush = () => {
		if (!block.some((line) => /@[0-9a-f]{6,}/.test(line))) out.push(...block);
		block = [];
	};
	for (const line of code.split('\n')) {
		if (line.startsWith('// ')) {
			// A printed value starts with `// name: Type =`; its continuation lines are indented.
			if (/^\/\/ [\w`$]+: /.test(line)) flush();
			block.push(line);
		} else {
			flush();
			if (!line.startsWith('assert(')) out.push(line);
		}
	}
	flush();
	return out.join('\n').replace(/\n{3,}/g, '\n\n').trim();
}

pages.forEach((page, order) => {
	const text = fs.readFileSync(path.join(compiled, page.source), 'utf8');
	const h1 = text.match(/^# (.+)$/m);
	if (!h1) fail(`docs/user/${page.source} has no H1 title`);
	const body = text.replace(h1[0], '').replace(/^\s+/, '');

	const rewritten = segments(body)
		.map((part) => (part.kind === 'code' ? tidyCode(part.text) : rewriteLinks(part.text, page.source)))
		.join('\n');

	// Anything still pointing at a relative .md file escaped the rewriter: a reference-style
	// link, a titled target, or raw HTML. Fail rather than publish a dead link.
	const leftover = rewritten.match(/\]\((?!https?:|\/|#)[^)]*\.md(?:#[^)]*)?\)|href="(?!https?:|\/|#)[^"]*\.md/);
	if (leftover) fail(`docs/user/${page.source} has a link the rewriter did not handle: ${leftover[0]}`);

	for (const snippet of landingSnippets.filter((s) => s.source === page.source)) {
		// An excerpt from the middle of a guide needs the guide's imports to stand alone.
		const imports = snippet.imports ? [...rewritten.matchAll(/^import .+$/gm)].map((m) => m[0]) : [];
		const blocks = snippet.headings.flatMap((heading) => {
			const section = rewritten.match(new RegExp(`^## ${heading}\\n([\\s\\S]*?)(?=^## |(?![\\s\\S]))`, 'm'));
			if (!section) fail(`docs/user/${page.source} no longer has a "${heading}" section; the landing page shows it`);
			return [...section[1].matchAll(/^```scala\n([\s\S]*?)\n```/gm)].map((m) => forLanding(m[1]));
		});
		snippets[snippet.name] = [...(imports.length > 0 ? [imports.join('\n')] : []), ...blocks].join('\n\n');
	}

	const description = describe(body);
	const frontmatter = [
		'---',
		`title: ${yaml(h1[1].trim())}`,
		description ? `description: ${yaml(description)}` : null,
		`editUrl: ${yaml(`${repo}/edit/${branch}/docs/user/${page.source}`)}`,
		'sidebar:',
		`  order: ${order}`,
		page.label ? `  label: ${yaml(page.label)}` : null,
		'---',
		'',
		// Markdown, not MDX: compiled Scala output is full of braces and angle brackets.
		`<!-- Generated by website/scripts/sync-docs.mjs from docs/user/${page.source}. Edit the source, not this file. -->`,
		'',
		'',
	].filter((line) => line !== null);

	const out = path.join(contentOut, `${page.slug}.md`);
	fs.mkdirSync(path.dirname(out), { recursive: true });
	fs.writeFileSync(out, frontmatter.join('\n') + rewritten);
});

fs.mkdirSync(generatedOut, { recursive: true });
fs.writeFileSync(path.join(generatedOut, 'snippets.json'), `${JSON.stringify(snippets, null, 2)}\n`);

console.log(`sync-docs: ${pages.length} pages from ${path.relative(repoRoot, compiled)}`);
