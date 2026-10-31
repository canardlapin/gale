#!/usr/bin/env node
// Draw the landing page's figure: the sparsity pattern of the five-point Laplacian on a
// k-by-k grid, beside the pattern of its Cholesky factor in the natural ordering. Entries
// of A are ink; entries that elimination creates (fill) are the accent colour.
//
// The figure is checked in as src/assets/cholesky-fill.svg. Rerun this script only to
// change it:   node scripts/draw-hero.mjs

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const k = 7;
const n = k * k;
const cell = 6;
const gap = 56;
const pad = 28;
const ink = '#1b2430';
const fill = '#2b7a8c';
const rule = '#c9cfd2';

// Pattern of A: each grid point couples to itself and its four neighbours.
const a = Array.from({ length: n }, () => new Array(n).fill(false));
for (let r = 0; r < k; r++) {
	for (let c = 0; c < k; c++) {
		const i = r * k + c;
		a[i][i] = true;
		if (c + 1 < k) a[i][i + 1] = a[i + 1][i] = true;
		if (r + 1 < k) a[i][i + k] = a[i + k][i] = true;
	}
}

// Symbolic Cholesky: eliminating column j joins every pair of rows below j that it touches.
const l = a.map((row, i) => row.map((v, j) => v && j <= i));
for (let j = 0; j < n; j++) {
	const below = [];
	for (let i = j + 1; i < n; i++) if (l[i][j]) below.push(i);
	for (const p of below) for (const q of below) if (q <= p) l[p][q] = true;
}

const nnzA = a.flat().filter(Boolean).length;
const nnzL = l.flat().filter(Boolean).length;
const nnzLowerA = a.flatMap((row, i) => row.filter((v, j) => v && j <= i)).length;

const side = n * cell;
const width = pad * 2 + side * 2 + gap;
const height = pad * 2 + side + 44;

const dots = (pattern, x0, colour) =>
	pattern
		.flatMap((row, i) =>
			row.map((v, j) => {
				if (!v) return '';
				const c = colour(i, j);
				return `<rect x="${x0 + j * cell + 0.6}" y="${pad + i * cell + 0.6}" width="${cell - 1.2}" height="${cell - 1.2}" rx="1" fill="${c}"/>`;
			}),
		)
		.join('');

const frame = (x0) => `<rect x="${x0 - 0.5}" y="${pad - 0.5}" width="${side + 1}" height="${side + 1}" fill="none" stroke="${rule}"/>`;
const xA = pad;
const xL = pad + side + gap;
const label = (x, text, sub) =>
	`<text x="${x + side / 2}" y="${pad + side + 34}" text-anchor="middle" font-family="Newsreader, Georgia, serif" font-size="24" fill="${ink}"><tspan font-style="italic">${text}</tspan><tspan font-family="Inter, system-ui, sans-serif" font-size="16" fill="#5b6670" dx="10">${sub}</tspan></text>`;

const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-labelledby="t d">
<title id="t">A sparse matrix and its Cholesky factor</title>
<desc id="d">Left: the ${nnzA} nonzeros of the five-point Laplacian on a ${k} by ${k} grid, a ${n} by ${n} matrix with five diagonals. Right: its Cholesky factor in the natural ordering, with ${nnzLowerA} entries inherited from A and ${nnzL - nnzLowerA} created by elimination, which fill the band.</desc>
${frame(xA)}${frame(xL)}
${dots(a, xA, () => ink)}
${dots(l, xL, (i, j) => (a[i][j] ? ink : fill))}
<text x="${xA + side + gap / 2}" y="${pad + side / 2 + 7}" text-anchor="middle" font-family="Inter, system-ui, sans-serif" font-size="20" fill="#5b6670">→</text>
${label(xA, 'A', `${nnzA} nonzeros`)}
${label(xL, 'L', `${nnzL - nnzLowerA} fill`)}
</svg>
`;

const out = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../src/assets/cholesky-fill.svg');
fs.writeFileSync(out, svg);
console.log(`draw-hero: n=${n}, nnz(A)=${nnzA}, nnz(L)=${nnzL}, fill=${nnzL - nnzLowerA}`);
