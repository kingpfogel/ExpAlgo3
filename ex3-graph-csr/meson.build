project('expalgo-csr', 'cpp')

omp = dependency('openmp')

executable('csr',
	'main.cpp',
	override_options: ['cpp_std=c++17'],
	dependencies: omp,
	install: true)
