// The site's page map. Every page comes from the mdoc output of docs/user, so the
// site publishes exactly what `sbt docs/mdoc` compiled and ran. Documents elsewhere
// under docs/ (plans, audits, release evidence) are not published here; links to
// them are rewritten to their GitHub URLs.

export const repo = 'https://github.com/canardlapin/gale';
export const branch = 'main';

/** @typedef {{ source: string, slug: string, label?: string }} Page */
/** @typedef {{ label: string, pages?: Page[], groups?: Group[], collapsed?: boolean }} Group */

/** @type {Group[]} */
export const groups = [
	{
		label: 'Start here',
		pages: [
			{ source: 'getting-started.md', slug: 'start/first-result', label: 'Your first result' },
			{ source: 'core-concepts.md', slug: 'start/core-concepts', label: 'Core concepts' },
			{ source: 'index.md', slug: 'start/reading-path', label: 'A reading path' },
		],
	},
	{
		label: 'Guides',
		pages: [
			{ source: 'guides/index.md', slug: 'guides', label: 'Choose a guide' },
			{ source: 'guides/examples.md', slug: 'guides/examples', label: 'Find a worked example' },
		],
		groups: [
			{
				label: 'Linear algebra',
				pages: [
					{ source: 'guides/dense-systems.md', slug: 'guides/dense-systems', label: 'Dense systems and least squares' },
					{ source: 'guides/reductions-and-numerics.md', slug: 'guides/reductions-and-numerics', label: 'Reductions and norms' },
					{ source: 'guides/centered-moments.md', slug: 'guides/centered-moments', label: 'Centered moments' },
					{ source: 'guides/sparse-operators.md', slug: 'guides/sparse-operators', label: 'Sparse matrices and operators' },
					{ source: 'guides/spectral-analysis.md', slug: 'guides/spectral-analysis', label: 'Spectral analysis' },
					{ source: 'guides/subspaces-and-minimum-norm.md', slug: 'guides/subspaces-and-minimum-norm', label: 'Numerical rank and subspaces' },
					{ source: 'guides/numerical-workspaces.md', slug: 'guides/compare-fits', label: 'Compare two linear fits' },
					{ source: 'guides/complex-hermitian.md', slug: 'guides/complex-hermitian', label: 'Complex Hermitian factors' },
				],
			},
			{
				label: 'Optimization',
				pages: [
					{ source: 'guides/first-order-optimization.md', slug: 'guides/first-order-optimization', label: 'First-order composite methods' },
					{ source: 'guides/nonlinear-and-bounded-optimization.md', slug: 'guides/nonlinear-and-bounded', label: 'Nonlinear fits and bounds' },
					{ source: 'guides/augmented-lagrangian.md', slug: 'guides/augmented-lagrangian', label: 'General constraints' },
					{ source: 'guides/cma-es.md', slug: 'guides/cma-es', label: 'Derivative-free (CMA-ES)' },
					{ source: 'guides/ica-rotation.md', slug: 'guides/ica-rotation', label: 'Symmetric tanh rotation' },
				],
			},
			{
				label: 'Migration',
				pages: [{ source: 'guides/breeze-equivalence.md', slug: 'guides/from-breeze', label: 'Moving from Breeze' }],
			},
			{
				label: 'Development previews',
				collapsed: true,
				pages: [{ source: 'guides/axis-margins.md', slug: 'guides/axis-margins', label: 'Axis margins' }],
			},
		],
	},
	{
		label: 'Advanced',
		pages: [
			{ source: 'advanced/index.md', slug: 'advanced', label: 'Overview' },
			{ source: 'advanced/numerical-contract.md', slug: 'advanced/numerical-contract', label: 'The numerical contract' },
			{ source: 'advanced/shipped-vs-deferred.md', slug: 'advanced/shipped-vs-deferred', label: 'Shipped and deferred' },
			{ source: 'advanced/generalized-operator-eigen.md', slug: 'advanced/generalized-eigen', label: 'Generalized eigensolving' },
			{ source: 'advanced/dense-destinations.md', slug: 'advanced/dense-destinations', label: 'Dense destinations' },
			{ source: 'advanced/workspaces.md', slug: 'advanced/workspaces', label: 'Reusable workspaces' },
			{ source: 'advanced/sparse-patterns.md', slug: 'advanced/sparse-patterns', label: 'Sparse patterns' },
			{ source: 'advanced/sparse-plans.md', slug: 'advanced/sparse-plans', label: 'Symbolic sparse plans' },
			{ source: 'advanced/ownership.md', slug: 'advanced/ownership', label: 'Ownership' },
		],
	},
	{
		label: 'Reference',
		pages: [
			{ source: 'reference/index.md', slug: 'reference', label: 'Overview' },
			{ source: 'reference/modules-and-platforms.md', slug: 'reference/modules-and-platforms', label: 'Modules and platforms' },
			{ source: 'reference/factorization-capabilities.md', slug: 'reference/factorization-capabilities', label: 'Factorization capabilities' },
			{ source: 'troubleshooting.md', slug: 'reference/troubleshooting', label: 'Troubleshooting' },
		],
	},
];

/** Every page in sidebar order. */
export function allPages() {
	/** @param {Group} g @returns {Page[]} */
	const walk = (g) => [...(g.pages ?? []), ...(g.groups ?? []).flatMap(walk)];
	return groups.flatMap(walk);
}

/** The Starlight sidebar for this map. */
export function sidebar() {
	/** @param {Group} g @returns {any} */
	const toItem = (g) => ({
		label: g.label,
		collapsed: g.collapsed ?? false,
		items: [
			...(g.pages ?? []).map((p) => (p.label ? { slug: p.slug, label: p.label } : p.slug)),
			...(g.groups ?? []).map(toItem),
		],
	});
	return groups.map(toItem);
}
